# For the stochastic volatility model with AR(1) dynamics (SVAR(1)), investigate the performance of efficient method of moments (EMM), rejection ABC and regression-IS ABC.

library(MASS)
library(mvtnorm)
library(tseries)
library(numDeriv)

# Simulation of the SVAR(1) model
randomness_sim<-function(sim_length){
	normal_error<-rnorm(sim_length)
	chisq_error<-rchisq(sim_length,1)
	signs<-sample(c(1,-1),sim_length,replace=T)
	randomerror<-cbind(normal_error,chisq_error,signs)
	return(randomerror)
}
SV_sim<-function(SVpars_randomerror,gen_random=FALSE){ # The input parameter is (\alpha,\beta,\sigma_{eta}), and the randomerror contains the normal errors, chisq errors and signs for y stored in columns; If N parameters are input in a matrix, the dimension should be 3*N; The output data is in a n*N matrix. 
	if(gen_random){
		randomerror<-randomness_sim(simu_length)
		SV_pars<-SVpars_randomerror
	}
	else{
		randomerror<-SVpars_randomerror[[2]]
		SV_pars<-SVpars_randomerror[[1]]
	}
	# Simulate through the lognormal SV model
	sim_vec<-c(SV_pars[2],abs(SV_pars[3]),SV_pars[1]/2/(1-SV_pars[2]))
	beta<-sim_vec[1]; sigma_eta<-sim_vec[2]; sigma_bar_star<-sim_vec[3]

	eta_t<-sigma_eta*randomerror[,1]
	x_t<-filter(eta_t,beta,method="recursive") # Simulate the AR(1) series
	epsilon_t_star<-randomerror[,2]
	y_t_star<-c(x_t)+log(epsilon_t_star)+2*sigma_bar_star
	y_t<-randomerror[,3]*exp(y_t_star/2)
	return(y_t)
}

# The MLE for the parameter of GARCH
GARCH_MLE_cov<-function(data){
	T_obs<-length(data)
	data_ts<-ts(data)
	garch_fit<-garch(data_ts,control=garch.control(grad="analytical",trace=F))
	garch_mle<-garch_fit$coef
	asycov<-garch_fit$vcov
	return(list(MLE=garch_mle,asycov=asycov))
}
GARCH_score<-function(data,garch_pars){
	# Calculate the score function of GARCH likelihood given the data and the parameters. First, the conditional variances are calculated; second, the gradient of the conditional variances are computed recursively; finally the gradient of the likelihood is calculated.

	a0<-garch_pars[1]; a1<-garch_pars[2]; b1<-garch_pars[3]
	# Compute the condition variance using the 'garch' function in 'tseries'
	T_data<-length(data)
	data_ts<-ts(data)
	garch_fit<-garch(data_ts,control=garch.control(grad="analytical",maxiter=1,start=garch_pars,trace=F))
	con_var<-(garch_fit$fitted.values^2)[,1]
	con_var[1]<-var(data)
	# Compute the gradient of the conditional variances
	grad_convar<-matrix(0,ncol=3,nrow=T_data)
	for(t_i in 2:T_data) grad_convar[t_i,]<-c(1,data[t_i-1]^2,con_var[t_i-1])+b1*grad_convar[t_i-1,]
	# Compute the score function
	score<-colSums(as.vector(data^2/con_var^2-1/con_var)*grad_convar)/2/T_data
	return(score)
}

# Transformation of the Simulation Parameters for Optimization
# Input and output a N*p matrix
transform_simupars<-function(simu_pars){ 
	if(is.vector(simu_pars)) return(c(simu_pars[1],log((1+simu_pars[2])/(1-simu_pars[2])),simu_pars[3]))
	else return(cbind(simu_pars[,1],log((1+simu_pars[,2])/(1-simu_pars[,2])),simu_pars[,3]))
}
invtransform_optimpars<-function(optimpars){
	if(is.vector(optimpars)) return(c(optimpars[1],(exp(optimpars[2])-1)/(exp(optimpars[2])+1),optimpars[3]))
	else return(cbind(optimpars[,1],(exp(optimpars[,2])-1)/(exp(optimpars[,2])+1),optimpars[,3]))
}

# Prior Density
priordens<-function(theta){
	N<-dim(theta)[1]; par_dim<-dim(theta)[2]
	prior_dens<-rep(1,N)
	return(prior_dens)
}
# Non-parametric Proposal Distribution
NP_prop<-function(cores,core_bandwidth,N){ # Output sample is a N*num_pars matrix
	N0<-dim(cores)[2]
	MC_mixture_ind<-sample(1:N0,N,replace=T)
	MC_mixture_error<-matrix(rnorm(num_pars*N,sd=core_bandwidth/sqrt(2*pi)),ncol=N) # num_pars*N matrix
	MC_NPprop<-MC_mixture_error+cores[,MC_mixture_ind]
	ttl_pairwise_sqdiff<-0
	for(p_i in 1:num_pars){
		pairwise_diff<-outer(MC_NPprop[p_i,],cores[p_i,],'-')/core_bandwidth[p_i] # N*N0 matrix
		ttl_pairwise_sqdiff<-ttl_pairwise_sqdiff+pairwise_diff^2
	}
	dens_NPprop<-exp(-pi*ttl_pairwise_sqdiff)%*%rep(1/N0,N0)
	MC_NPprop<-t(MC_NPprop) # N*num_pars matrix
	return(list(sample=MC_NPprop,density_prop=dens_NPprop))
}

###################################################
# Summary set 1 (Using Efficient Method of Moments)
###################################################

# The Efficient Method of Moments
EMM_crit<-function(optim_pars,randomerror,obs_auxpars,inv_weight_matrix){
	SV_pars<-invtransform_optimpars(optim_pars)
	if(abs(SV_pars[2])>.99) return(Inf)
	simu_data<-SV_sim(list(SV_pars,randomerror))
	simu_score<-GARCH_score(simu_data,obs_auxpars)
	criteria<-c(t(simu_score)%*%inv_weight_matrix%*%simu_score)
	return(criteria)
}
EMM<-function(randomerror,obs_auxpars,inv_weight_matrix){
	optim_parsini<-transform_simupars(true_pars)
	EMM_results<-optim(optim_parsini,EMM_crit,method='BFGS',randomerror=randomerror,obs_auxpars=obs_auxpars,inv_weight_matrix=inv_weight_matrix,control=list(trace=T))
	EMM_pars<-invtransform_optimpars(EMM_results$par)
	names(EMM_pars)<-list('alpha','beta','sigma')
	return(EMM_pars)
}

# Summary Statistics
ABC_summary<-function(SVpars,gen_random=TRUE,randomerror=NULL){
	if(gen_random){
		simu_data<-SV_sim(SVpars,gen_random=TRUE)
		simu_results<-GARCH_MLE_cov(simu_data)
		sum_stat<-simu_results$MLE
	}
	if(!gen_random){
		SVpars_randomerror<-list(SVpars,randomerror)
		simu_data<-SV_sim(SVpars_randomerror,gen_random=FALSE)
		simu_results<-GARCH_MLE_cov(simu_data)
		sum_stat<-simu_results$MLE
	}
	return(sum_stat)
}

# Regression IS-ABC
RIS_ABC<-function(observations,cores,cores_trans_bandwidth,N){ # core is a N0*num_pars matrix
	cores_trans<-transform_simupars(cores)
	MC_results<-NP_prop(t(cores_trans),cores_trans_bandwidth,N)
	MC_sample<-MC_results$sample # N*num_pars matrix
	MC_propdens<-MC_results$density_prop
	MC_pridens<-priordens(MC_sample)
	MC_sample_invtrans<-invtransform_optimpars(MC_sample)
	obs_summary<-GARCH_MLE_cov(observations)$MLE
	MC_summmary<-apply(MC_sample_invtrans,1,ABC_summary) # num_auxpars*N matrix
	MC_diff<-MC_summmary-obs_summary
	MC_diff<-t(MC_diff) # N*num_auxpars matrix 

	MC_finind<-is.finite(rowSums(MC_diff))
	MC_diff<-MC_diff[MC_finind,]
	MC_sample_invtrans<-MC_sample_invtrans[MC_finind,]
	MC_propdens<-MC_propdens[MC_finind]
	MC_pridens<-MC_pridens[MC_finind]

	MC_regresults<-lm(MC_sample_invtrans~MC_diff-1)
	MC_regcoef<-MC_regresults$coef
	MC_linadjused_sample<-MC_sample_invtrans-MC_diff%*%MC_regcoef
    MC_ISweights<-c(MC_pridens/MC_propdens)
    ABC_est<-colSums(MC_linadjused_sample*MC_ISweights)/sum(MC_ISweights)
    return(list(est=ABC_est,sample=MC_sample,diff=MC_diff,weights=MC_ISweights))
}

# Setup
set.seed(100)
true_pars<-c(-0.736,0.9,0.363); T_obs<-500; num_pars<-length(true_pars); num_auxpars<-3
obs_randomerror<-cbind(rnorm(T_obs),rchisq(T_obs,1),sample(c(1,-1),T_obs,replace=T))
observations<-SV_sim(list(true_pars,obs_randomerror))

prior_alpha_lower<--4; prior_alpha_upper<-2
prior_beta_lower<--2; prior_beta_upper<-2
prior_sigmaeta_lower<--1; prior_sigmaeta_upper<-2
prior_unifpar<-matrix(c(prior_alpha_lower,prior_beta_lower,prior_sigmaeta_lower,prior_alpha_upper,prior_beta_upper,prior_sigmaeta_upper),num_pars,2)


# MC distribution of EMM estimates and the estimated posterior distribution
obs_auxresults<-GARCH_MLE_cov(observations)
obs_auxpars<-obs_auxresults$MLE
weight_matrix<-solve(obs_auxresults$asycov*T_obs)
inv_weight_matrix<-solve(weight_matrix)
MC_EMM_detjacob<-0
simu_length<-20000; N0<-1000; MC_EMM<-matrix(0,N0,num_pars)

x11()
plot(0,0,xlim=c(1,N0),ylim=c(1,N0))
set.seed(100)
for(N_i in 1:N0){
	points(N_i,N_i,pch=20)
	simu_randomerror<-randomness_sim(sim_length=simu_length)
	options(show.error.messages = FALSE)
	system.time(temp<-try(EMM(simu_randomerror,obs_auxpars,inv_weight_matrix)))
	if(!is.numeric(temp[1])) {N_i<-N_i-1;next}
	MC_EMM[N_i,]<-temp
	MC_EMM_detjacob[N_i]<-abs(det(jacobian(ABC_summary,MC_EMM[N_i,],gen_random=FALSE,randomerror=simu_randomerror)))
}
options(show.error.messages = TRUE)

# Reverse sampler of Forneron and Ng(2015) and Meeds and Welling(2015)
finite.ind<-is.finite(MC_EMM_detjacob)
MC_EMM_detjacob<-MC_EMM_detjacob[finite.ind]
MC_EMM<-MC_EMM[finite.ind,]
N0<-sum(finite.ind)
EMM_pridens<-priordens(MC_EMM)
EMM_ISweights<-c(EMM_pridens/MC_EMM_detjacob)
EMM_postmean<-colSums(MC_EMM*EMM_ISweights)/sum(EMM_ISweights)

replic<-floor(N0/50); EMM_postmean_all<-matrix(0,replic,num_pars); ind_all<-1:N0
for(rep_i in 1:replic){
	rep_i_ind<-sample(ind_all,50)#((rep_i-1)*50+1):(rep_i*50)
	ind_all<-setdiff(ind_all,rep_i_ind)
	EMM_postmean_all[rep_i,]<-colSums(MC_EMM[rep_i_ind,]*EMM_ISweights[rep_i_ind])/sum(EMM_ISweights[rep_i_ind])
}
EMM_postmean_MSE<-rowMeans((t(EMM_postmean_all)-true_pars)^2)
EMM_postmean_bias2<-(colMeans(EMM_postmean_all)-true_pars)^2
EMM_postmean_var<-apply(EMM_postmean_all,2,var)
signif(cbind(MSE=EMM_postmean_MSE,bias2=EMM_postmean_bias2,var=EMM_postmean_var),2)

# Effective sample size of the importance sampling weights of the reverse sampler
rep_ind<-sample(1:N0,N0)
1/(1+var(EMM_ISweights[rep_ind])/mean(EMM_ISweights[rep_ind])^2)
weight_order<-order(EMM_ISweights,decreasing=T)
EMM_ISweights[weight_order[1:10]]
MC_EMM[weight_order[1:10],]

# Kernel density estimates of the proposal density and the posterior density
xinterv<-cbind(c(-4,0),c(0.5,1),c(-1,1.5))
yinterv<-cbind(c(0,4),c(0,30),c(0,12))
ind_all<-1:N0
rep_ind<-sample(ind_all,N0)
x11()
par(mfrow=c(2,2))
for(plot_i in 1:num_pars){
	EMM_dens<-density(MC_EMM[rep_ind,plot_i])
	EMM_post_dens<-density(MC_EMM[rep_ind,plot_i],weights=EMM_ISweights[rep_ind]/sum(EMM_ISweights[rep_ind]))
	ymax<-max(EMM_dens$y,EMM_post_dens$y)
	ymin<-min(EMM_dens$y,EMM_post_dens$y)
	plot(EMM_dens,ylim=yinterv[,plot_i],xlim=xinterv[,plot_i],main='Reverse Sampler summary set 1')
	lines(EMM_post_dens,ylim=c(ymin,ymax),col=3)
	abline(v=EMM_postmean[plot_i],col=3,lty=2)
	abline(v=mean(MC_EMM[,plot_i]),col=1,lty=2)
	abline(v=true_pars[plot_i],col=2,lty=2)
}

# RegIS-ABC
obs_auxresults<-GARCH_MLE_cov(observations)
obs_auxpars<-obs_auxresults$MLE
weight_matrix<-solve(obs_auxresults$asycov*T_obs)
inv_weight_matrix<-solve(weight_matrix)
simu_length<-500; N0<-50; MC_EMM<-matrix(0,N0,num_pars)

set.seed(100)
for(N_i in 1:N0){
	simu_randomerror<-randomness_sim(sim_length=simu_length)
	MC_EMM[N_i,]<-EMM(simu_randomerror,obs_auxpars,inv_weight_matrix)
}

N<-300; replic<-300; 
RegIS_ABC_all<-matrix(0,replic,num_pars)
RegIS_ABC_sample<-list(0)
RegIS_ABC_diff<-list(0)
RegIS_ABC_weight<-list(0)

set.seed(100)
MC_EMM_trans<-transform_simupars(MC_EMM)
MC_EMM_trans_bandwidth<-apply(MC_EMM_trans,2,sd)/10#*sqrt(2*pi)
for(rep_i in 1:replic){
	RegIS_ABC_results<-RIS_ABC(observations,MC_EMM,MC_EMM_trans_bandwidth,N)
	RegIS_ABC_all[rep_i,]<-RegIS_ABC_results$est
	RegIS_ABC_sample[[rep_i]]<-RegIS_ABC_results$sample
	RegIS_ABC_diff[[rep_i]]<-RegIS_ABC_results$diff
	RegIS_ABC_weight[[rep_i]]<-RegIS_ABC_results$weights
}

acc_quant<-2/3
ABC1_est<-matrix(0,replic,num_pars); ABC2_est<-matrix(0,replic,num_pars); ABC3_est<-matrix(0,replic,num_pars); ABC4_est<-matrix(0,replic,num_pars); ABC5_est<-matrix(0,replic,num_pars)
for(rep_i in 1:replic){
	MC_sample<-RegIS_ABC_sample[[rep_i]]
	MC_diff<-RegIS_ABC_diff[[rep_i]]
	MC_dist<-sqrt(rowSums(MC_diff^2))
    MC_ISweights<-RegIS_ABC_weight[[rep_i]]
	MC_accind<-(MC_dist<quantile(MC_dist,acc_quant))
	MC_sample<-MC_sample[MC_accind,]
	MC_diff<-MC_diff[MC_accind,]
	MC_ISweights<-MC_ISweights[MC_accind]

	# Regression of parameter values over distances
	MC_regresults1<-lm(MC_sample~MC_diff)
	MC_regcoef1<-MC_regresults1$coef
	MC_linadjused_sample1<-MC_sample-MC_diff%*%MC_regcoef1[-1,]
	MC_linadjused_sample2<-cbind(1,MC_diff)%*%MC_regcoef1
	MC_diff2<-poly(MC_diff,degree=2,raw=T)
	MC_regresults2<-lm(MC_sample~poly(MC_diff,degree=2,raw=T))
	MC_regcoef2<-MC_regresults2$coef
	MC_quadradjused_sample1<-MC_sample-MC_diff2%*%MC_regcoef2[-1,]
	MC_quadradjused_sample2<-cbind(1,MC_diff2)%*%MC_regcoef2

    ABC1_est[rep_i,]<-colSums(invtransform_optimpars(MC_linadjused_sample1)*MC_ISweights)/sum(MC_ISweights) # Estimator using linear predictions
    ABC2_est[rep_i,]<-colSums(invtransform_optimpars(MC_linadjused_sample2)*MC_ISweights)/sum(MC_ISweights) # Estimator using linear condition expectations
    ABC3_est[rep_i,]<-colSums(invtransform_optimpars(MC_sample)*MC_ISweights)/sum(MC_ISweights) # No regression
    ABC4_est[rep_i,]<-colSums(invtransform_optimpars(MC_quadradjused_sample1)*MC_ISweights)/sum(MC_ISweights) # Estimator using quadratic predictions
    ABC5_est[rep_i,]<-colSums(invtransform_optimpars(MC_quadradjused_sample2)*MC_ISweights)/sum(MC_ISweights) # Estimator using quadratic condition expectations
}
ABC1_MSE<-rowMeans((t(ABC1_est)-true_pars)^2)
ABC1_bias2<-(colMeans(ABC1_est)-true_pars)^2
ABC1_var<-apply(ABC1_est,2,var)
ABC2_MSE<-rowMeans((t(ABC2_est)-true_pars)^2)
ABC2_bias2<-(colMeans(ABC2_est)-true_pars)^2
ABC2_var<-apply(ABC2_est,2,var)
ABC3_MSE<-r 4_MSE<-rowMeans((t(ABC4_est)-true_pars)^2)
ABC4_bias2<-(colMeans(ABC4_est)-true_pars)^2
ABC4_var<-apply(ABC4_est,2,var)
ABC5_MSE<-rowMeans((t(ABC5_est)-true_pars)^2)
ABC5_bias2<-(colMeans(ABC5_est)-true_pars)^2
ABC5_var<-apply(ABC5_est,2,var)

signif(cbind(ABC1_MSE,ABC2_MSE,ABC3_MSE,ABC4_MSE,ABC5_MSE),2)
signif(cbind(ABC1_bias2,ABC2_bias2,ABC3_bias2,ABC4_bias2,ABC5_bias2),2)
signif(cbind(ABC1_var,ABC2_var,ABC3_var,ABC4_var,ABC5_var),2)


rep_i<-1; set.seed(103)
cores_trans<-transform_simupars(cores)
cores_trans_bandwidth<-apply(cores_trans,2,sd)/10#*sqrt(2*pi)
RegIS_ABC_results<-RIS_ABC(observations,MC_EMM,N)
RegIS_ABC_all[rep_i,]<-RegIS_ABC_results$est
RegIS_ABC_sample[[rep_i]]<-RegIS_ABC_results$sample
RegIS_ABC_diff[[rep_i]]<-RegIS_ABC_results$diff
RegIS_ABC_weight[[rep_i]]<-RegIS_ABC_results$weights

acc_quant<-2/3
MC_sample<-RegIS_ABC_sample[[rep_i]]
MC_diff<-RegIS_ABC_diff[[rep_i]]
MC_dist<-sqrt(rowSums(MC_diff^2))
MC_ISweights<-RegIS_ABC_weight[[rep_i]]
MC_accind<-(MC_dist<quantile(MC_dist,acc_quant))
MC_sample<-MC_sample[MC_accind,]
MC_diff<-MC_diff[MC_accind,]
MC_ISweights<-MC_ISweights[MC_accind]
# MC_regresults<-lm(MC_sample~cbind(MC_diff,MC_diff^2))
MC_regresults<-lm(MC_sample~MC_diff)
MC_regcoef<-MC_regresults$coef
MC_linadjused_sample1<-MC_sample-MC_diff%*%MC_regcoef[-1,]
# MC_linadjused_sample<-cbind(1,MC_diff,MC_diff^2)%*%MC_regcoef
MC_linadjused_sample<-cbind(1,MC_diff)%*%MC_regcoef
ABC<-colSums(invtransform_optimpars(MC_sample)*MC_ISweights)/sum(MC_ISweights)
ABC_linadjest<-colSums(invtransform_optimpars(MC_linadjused_sample)*MC_ISweights)/sum(MC_ISweights)

xinterv<-cbind(c(-4,0),c(0.5,1),c(-1,1.5))
yinterv<-cbind(c(0,4),c(0,30),c(0,12))
x11()
par(mfrow=c(2,2))
for(plot_i in 1:num_pars){
	EMM_post_dens<-density(MC_EMM[,plot_i],weights=EMM_ISweights/sum(EMM_ISweights))
	MC_linadjust_dens<-density(invtransform_optimpars(MC_linadjused_sample)[,plot_i])
	MC_linadjust_dens1<-density(invtransform_optimpars(MC_linadjused_sample1)[,plot_i])
	MC_dens<-density(invtransform_optimpars(MC_sample)[,plot_i])	

	MC_linadjust_post_dens<-density(invtransform_optimpars(MC_linadjused_sample)[,plot_i],weights=MC_ISweights/sum(MC_ISweights))
	MC_linadjust_post_dens1<-density(invtransform_optimpars(MC_linadjused_sample1)[,plot_i],weights=MC_ISweights/sum(MC_ISweights))
	MC_post_dens<-density(invtransform_optimpars(MC_sample)[,plot_i],weights=MC_ISweights/sum(MC_ISweights))

	xmax<-max(c(MC_EMM[,plot_i],invtransform_optimpars(MC_sample)[,plot_i],invtransform_optimpars(MC_linadjused_sample)[,plot_i]))
	xmin<-min(c(MC_EMM[,plot_i],invtransform_optimpars(MC_sample)[,plot_i],invtransform_optimpars(MC_linadjused_sample)[,plot_i]))
	ymax<-max(EMM_post_dens$y,MC_linadjust_dens$y,MC_post_dens$y,MC_dens$y)
	ymin<-min(EMM_post_dens$y,MC_linadjust_dens$y,MC_post_dens$y,MC_dens$y)
	plot(EMM_post_dens,xlim=xinterv[,plot_i],ylim=yinterv[,plot_i],main=paste('Accept',signif(acc_quant,2)))
	lines(MC_linadjust_dens,col=2)
	lines(MC_linadjust_dens1,col=2,lty=2)
	lines(MC_dens,col=3)
	# abline(v=EMM_postmean[plot_i],col=1,lty=2)	
	abline(v=true_pars[plot_i],col=2,lty=2)
	# abline(v=ABC_linadjest[plot_i],col=3,lty=2)
	# abline(v=ABC[plot_i],col=4,lty=2)	
}



xinterv<-cbind(c(-4,0),c(0.5,1),c(-1,1.5))
yinterv<-cbind(c(0,4),c(0,30),c(0,12))
x11()
par(mfrow=c(2,2))
for(plot_i in 1:num_pars){
	EMM_dens<-density(MC_EMM[,plot_i])
	EMM_post_dens<-density(MC_EMM[,plot_i],weights=EMM_ISweights/sum(EMM_ISweights))
	MC_dens<-density(invtransform_optimpars(MC_sample)[,plot_i])	
	MC_post_dens<-density(invtransform_optimpars(MC_sample)[,plot_i],weights=MC_ISweights/sum(MC_ISweights))
	ymax<-max(EMM_dens$y,MC_post_dens$y,MC_dens$y,EMM_post_dens$y)
	ymin<-min(EMM_dens$y,MC_post_dens$y,MC_dens$y,EMM_post_dens$y)
	plot(EMM_dens,xlim=xinterv[,plot_i],ylim=yinterv[,plot_i],main=paste('Accept',signif(acc_quant,2)))
	lines(EMM_post_dens,lty=2)
	lines(MC_dens,col=4)
	lines(MC_post_dens,col=4,lty=2)
	abline(v=true_pars[plot_i],col=2,lty=2)
}


###################################################
# Summary set 2 (Using Indirect Inference)
###################################################

# In the first stage, instead of using the EMM estimates, use the parameter values obtained from optimizing the following summary statistic.
cal_summary<-function(data){
	n<-length(data)
	s3<-t(rep(1/n,n))%*%data
	s2<-t(rep(1/n,n))%*%data^2-s3^2
	s1<-(t(rep(1/(n-1),n-1))%*%(data[-1]*data[-n])-s3^2)/s2
	return(c(s1,s2,s3))
}
ABC_summary<-function(SVpars,gen_random=TRUE,randomerror=NULL){
	if(gen_random){
		simu_data<-SV_sim(SVpars,gen_random=TRUE)
		simu_sum<-cal_summary(simu_data)
	}
	if(!gen_random){
		SVpars_randomerror<-list(SVpars,randomerror)
		simu_data<-SV_sim(SVpars_randomerror,gen_random=FALSE)
		simu_sum<-cal_summary(simu_data)
	}
	return(simu_sum)
}

sum_criteria<-function(optim_pars,randomerror,obs_sum){
	SV_pars<-invtransform_optimpars(optim_pars)
	if(abs(SV_pars[2])>.99||!is.finite(SV_pars[2])) return(Inf)
	simu_sum<-ABC_summary(SV_pars,gen_random=FALSE,randomerror=randomerror)
	criteria<-sqrt(sum((simu_sum-obs_sum)^2))
	return(criteria)
}
sum_optim<-function(randomerror,observations){
	obs_sum<-cal_summary(observations)
	optim_parsini<-transform_simupars(true_pars)
	optim_results<-optim(optim_parsini,sum_criteria,method='BFGS',randomerror=randomerror,obs_sum=obs_sum,control=list(trace=T))
	optim_pars<-invtransform_optimpars(optim_results$par)
	names(optim_pars)<-list('alpha','beta','sigma')
	return(optim_pars)
}


# Regression IS-ABC
RIS_ABC<-function(observations,cores,cores_trans_bandwidth,N){ # core is a N0*num_pars matrix
	cores_trans<-transform_simupars(cores)
	MC_results<-NP_prop(t(cores_trans),cores_trans_bandwidth,N)
	MC_sample<-MC_results$sample # N*num_pars matrix
	MC_propdens<-MC_results$density_prop
	MC_pridens<-priordens(MC_sample)
	MC_sample_invtrans<-invtransform_optimpars(MC_sample)
	obs_summary<-cal_summary(observations)
	MC_summmary<-apply(MC_sample_invtrans,1,ABC_summary) # num_auxpars*N matrix
	MC_diff<-MC_summmary-obs_summary
	MC_diff<-t(MC_diff) # N*num_auxpars matrix 

	MC_finind<-is.finite(rowSums(MC_diff))
	MC_diff<-MC_diff[MC_finind,]
	MC_sample_invtrans<-MC_sample_invtrans[MC_finind,]
	MC_propdens<-MC_propdens[MC_finind]
	MC_pridens<-MC_pridens[MC_finind]

	MC_regresults<-lm(MC_sample_invtrans~MC_diff-1)
	MC_regcoef<-MC_regresults$coef
	MC_linadjused_sample<-MC_sample_invtrans-MC_diff%*%MC_regcoef
    MC_ISweights<-c(MC_pridens/MC_propdens)
    ABC_est<-colSums(MC_linadjused_sample*MC_ISweights)/sum(MC_ISweights)
    return(list(est=ABC_est,sample=MC_sample,diff=MC_diff,weights=MC_ISweights))
}

# Setup
set.seed(100)
true_pars<-c(-0.736,0.9,0.363); T_obs<-2000; num_pars<-length(true_pars); num_auxpars<-3
obs_randomerror<-cbind(rnorm(T_obs),rchisq(T_obs,1),sample(c(1,-1),T_obs,replace=T))
observations<-SV_sim(list(true_pars,obs_randomerror))

prior_alpha_lower<--4; prior_alpha_upper<-2
prior_beta_lower<--2; prior_beta_upper<-2
prior_sigmaeta_lower<--1; prior_sigmaeta_upper<-2
prior_unifpar<-matrix(c(prior_alpha_lower,prior_beta_lower,prior_sigmaeta_lower,prior_alpha_upper,prior_beta_upper,prior_sigmaeta_upper),num_pars,2)


# MC distribution of II estimates and the estimated posterior distribution
simu_length<-T_obs; N0<-50
MC_optimsum<-matrix(0,N0,num_pars); MC_optimsum_detjacob<-0
set.seed(100)
for(N_i in 1:N0){
	simu_randomerror<-randomness_sim(sim_length=simu_length)
	options(show.error.messages = FALSE)
	temp<-try(sum_optim(simu_randomerror,observations))
	if(!is.numeric(temp[1])) {N_i<-N_i-1;next}
	MC_optimsum[N_i,]<-temp
	MC_optimsum_detjacob[N_i]<-abs(det(jacobian(ABC_summary,MC_optimsum[N_i,],gen_random=FALSE,randomerror=simu_randomerror)))
}
options(show.error.messages = TRUE)

# Reverse sampler of Forneron and Ng(2015) and Meeds and Welling(2015)
finite.ind<-is.finite(MC_optimsum_detjacob)
MC_optimsum_detjacob<-MC_optimsum_detjacob[finite.ind]
MC_optimsum<-MC_optimsum[finite.ind,]
N0<-sum(finite.ind)
optimsum_pridens<-priordens(MC_optimsum)
optimsum_ISweights<-c(optimsum_pridens/MC_optimsum_detjacob)
optimsum_postmean<-colSums(MC_optimsum*optimsum_ISweights)/sum(optimsum_ISweights)

replic<-floor(N0/50); optimsum_postmean_all<-matrix(0,replic,num_pars); ind_all<-1:N0
for(rep_i in 1:replic){
	rep_i_ind<-sample(ind_all,50)#((rep_i-1)*50+1):(rep_i*50)
	ind_all<-setdiff(ind_all,rep_i_ind)
	optimsum_postmean_all[rep_i,]<-colSums(MC_optimsum[rep_i_ind,]*optimsum_ISweights[rep_i_ind])/sum(optimsum_ISweights[rep_i_ind])
}
optimsum_postmean_MSE<-rowMeans((t(optimsum_postmean_all)-true_pars)^2)
optimsum_postmean_bias2<-(colMeans(optimsum_postmean_all)-true_pars)^2
optimsum_postmean_var<-apply(optimsum_postmean_all,2,var)
signif(cbind(MSE=optimsum_postmean_MSE,bias2=optimsum_postmean_bias2,var=optimsum_postmean_var),2)

# Effective sample size of the importance sampling weights of the reverse sampler
optimsum_pridens<-priordens(MC_optimsum)
optimsum_ISweights<-c(optimsum_pridens/MC_optimsum_detjacob)
rep_ind<-sample(1:N0,N0)
1/(1+var(optimsum_ISweights[rep_ind])/mean(optimsum_ISweights[rep_ind])^2)
weight_order<-order(optimsum_ISweights,decreasing=T)
optimsum_ISweights[weight_order[1:10]]
MC_optimsum[weight_order[1:10],]

# Kernel density estimates of the proposal density and the posterior density
xinterv<-cbind(c(-4,0),c(0.5,1),c(-1,1.5))
yinterv<-cbind(c(0,4),c(0,30),c(0,12))
ind_all<-1:N0
rep_ind<-sample(ind_all,50)
x11()
par(mfrow=c(2,2))
for(plot_i in 1:num_pars){
	optimsum_dens<-density(MC_optimsum[rep_ind,plot_i])
	optimsum_post_dens<-density(MC_optimsum[rep_ind,plot_i],weights=optimsum_ISweights[rep_ind]/sum(optimsum_ISweights[rep_ind]))
	ymax<-max(optimsum_dens$y,optimsum_post_dens$y)
	ymin<-min(optimsum_dens$y,optimsum_post_dens$y)
	plot(optimsum_dens,ylim=yinterv[,plot_i],xlim=xinterv[,plot_i])
	lines(optimsum_post_dens,ylim=c(ymin,ymax),col=3)
	abline(v=optimsum_postmean[plot_i],col=3,lty=2)
	abline(v=true_pars[plot_i],col=2,lty=2)
}

# RegIS-ABC
N<-300; replic<-100; 
RegIS_ABC_all<-matrix(0,replic,num_pars)
RegIS_ABC_sample<-list(0)
RegIS_ABC_diff<-list(0)
RegIS_ABC_weight<-list(0)

set.seed(100)
MC_optimsum_trans<-transform_simupars(MC_optimsum)
MC_optimsum_trans_bandwidth<-apply(MC_optimsum_trans,2,sd)/10#*sqrt(2*pi)
for(rep_i in 1:replic){
	RegIS_ABC_results<-RIS_ABC(observations,MC_optimsum,MC_optimsum_trans_bandwidth,N)
	RegIS_ABC_all[rep_i,]<-RegIS_ABC_results$est
	RegIS_ABC_sample[[rep_i]]<-RegIS_ABC_results$sample
	RegIS_ABC_diff[[rep_i]]<-RegIS_ABC_results$diff
	RegIS_ABC_weight[[rep_i]]<-RegIS_ABC_results$weights
}

rep_i<-3
1/(1+var(RegIS_ABC_weight[[rep_i]])/mean(RegIS_ABC_weight[[rep_i]])^2)

acc_quant<-1/10
ABC1_est<-matrix(0,replic,num_pars); ABC2_est<-matrix(0,replic,num_pars); ABC3_est<-matrix(0,replic,num_pars)
for(rep_i in 1:replic){
	MC_sample<-RegIS_ABC_sample[[rep_i]]
	MC_diff<-RegIS_ABC_diff[[rep_i]]
	MC_dist<-sqrt(rowSums(MC_diff^2))
    MC_ISweights<-RegIS_ABC_weight[[rep_i]]
	MC_accind<-(MC_dist<quantile(MC_dist,acc_quant))
	MC_sample<-MC_sample[MC_accind,]
	MC_diff<-MC_diff[MC_accind,]
	MC_ISweights<-MC_ISweights[MC_accind]

	# Regression of parameter values over distances
	MC_regresults<-lm(MC_sample~MC_diff)
	MC_regcoef<-MC_regresults$coef
	MC_linadjused_sample1<-MC_sample-MC_diff%*%MC_regcoef[-1,]
	MC_linadjused_sample2<-cbind(1,MC_diff)%*%MC_regcoef
    ABC1_est[rep_i,]<-colSums(invtransform_optimpars(MC_linadjused_sample1)*MC_ISweights)/sum(MC_ISweights) # Estimator using predictions
    ABC2_est[rep_i,]<-colSums(invtransform_optimpars(MC_linadjused_sample2)*MC_ISweights)/sum(MC_ISweights) # Estimator using condition expectations
    ABC3_est[rep_i,]<-colSums(invtransform_optimpars(MC_sample)*MC_ISweights)/sum(MC_ISweights) # No regression

}
ABC1_MSE<-rowMeans((t(ABC1_est)-true_pars)^2)
ABC1_bias2<-(colMeans(ABC1_est)-true_pars)^2
ABC1_var<-apply(ABC1_est,2,var)
ABC2_MSE<-rowMeans((t(ABC2_est)-true_pars)^2)
ABC2_bias2<-(colMeans(ABC2_est)-true_pars)^2
ABC2_var<-apply(ABC2_est,2,var)
ABC3_MSE<-rowMeans((t(ABC3_est)-true_pars)^2)
ABC3_bias2<-(colMeans(ABC3_est)-true_pars)^2
ABC3_var<-apply(ABC3_est,2,var)

signif(cbind(ABC1_MSE,ABC2_MSE,ABC3_MSE),2)
signif(cbind(ABC1_bias2,ABC2_bias2,ABC3_bias2),2)
signif(cbind(ABC1_var,ABC2_var,ABC3_var),2)

