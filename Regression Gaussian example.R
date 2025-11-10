library(mcmc)
load("/Volumes/Data_Drive/normquantile_pool2.RData")
load('/Volumes/Data_Drive/n1e4_experiment.RData')


######################################################################################################
# Pre-generate the pool of normal sample quantiles, since the maximum data size*MCsize=10^5*10^6 is too large 
# to re-simulate in every replicated experiment
######################################################################################################
# All quantiles are generated together and each tested quantile set is taken out when required
random_seed<-100; set.seed(random_seed)
n_all<-10^(2:5); n_max<-max(n_all); n_length<-length(n_all)
N_max<-10^6; N_max_part<-10^8/n_max; N_part_no<-N_max/N_max_part
probs_size_all<-c(2,4,9,19); probs_length<-length(probs_size_all)
# Index sets for each tested quantile set
probs_all_indices<-list(0)
for(i in 1:probs_length){
	if(i==1) probs_all_indices[[i+1]]<-(probs_all_indices[[i]][1]+1):(probs_all_indices[[i]][1]+probs_size_all[i])
	if(i>1) probs_all_indices[[i+1]]<-(probs_all_indices[[i]][probs_size_all[i-1]]+1):(probs_all_indices[[i]][probs_size_all[i-1]]+probs_size_all[i])
}
# Generate normal quantiles; output normquantile_all which is N_max*sum(probs_size_all) matrix
probs_all<-list(0)
for(probs_i in 1:probs_length) {
	probs_size<-probs_size_all[probs_i]; 
	probs_all[[probs_i]]<-seq(0,1,length=probs_size_all[probs_i]+2)[c(-1,-(probs_size+2))]
}
probs_all_unlist<-unlist(probs_all)
normquantile_all<-list(NULL,NULL,NULL,NULL)
x11(); plot(0,0,xlim=c(1,N_part_no),ylim=c(1,N_part_no))
for(i in 1:N_part_no){
	points(i,i,pch=20)
	normdata<-matrix(rnorm(n_max*N_max_part),nrow=N_max_part,ncol=n_max)
	for(n_i in 1:n_length){
		n<-n_all[n_i]
		normquantile<-t(apply(normdata[,1:n],1,quantile,probs=probs_all_unlist)) # sum(probs_size_all)*N_max_part
		normquantile_all[[n_i]]<-rbind(normquantile_all[[n_i]],normquantile)		
	} 
}
dev.off()
probs_indices_pool<-probs_indices[-1]
normquantile_pool<-normquantile_all
rm(list=setdiff(ls(), c("normquantile_pool","probs_indices_pool","probs_size_all")))
save.image("E:\\Dropbox\\Coding\\normquantile_pool-Gaussian.RData")



##############################################################################
# Functions
##############################################################################
rNormal<-function(n,mu,sigma){ # Simulating univariate normal observations
	observations<-rnorm(n,mu,sigma)
	return(observations)
}
negloglik<-function(theta,sumstat){ # Output the negative log-likelihood function of the summary statistics on the given parameter value theta
	mu<-theta[1]; sigma<-abs(theta[2])
	quant<-sign(sumstat)*sqrt(abs(sumstat)); d<-length(sumstat)
	# quant<-2*log(sumstat); d<-length(quant)
	# quant<-sumstat; d<-length(sumstat)
	quant_norm<-(quant-mu)/sigma
	pnorm_quant<-pnorm(quant_norm)
	pnorm_diff<-c(pnorm_quant,1)-c(0,pnorm_quant)
	n_interv_adj<-n_interv-1
	negloglik_val<-sum(quant_norm^2)/2+d*log(sigma)-sum(n_interv_adj*log(pnorm_diff))
	return(negloglik_val)
}
neglikscore<-function(theta,sumstat){ # Output the score function for the negative log-likelihood
	mu<-theta[1]; sigma<-abs(theta[2])
	quant<-sign(sumstat)*sqrt(abs(sumstat)); d<-length(sumstat)
	# quant<-2*log(sumstat); d<-length(sumstat)
	# quant<-sumstat; d<-length(sumstat)
	quant_norm<-(quant-mu)/sigma
	pnorm_quant<-pnorm(quant_norm)
	pnorm_diff<-c(pnorm_quant,1)-c(0,pnorm_quant)
	dnorm_quant<-dnorm(quant_norm)	
	dnorm_diff<-c(dnorm_quant,0)-c(0,dnorm_quant)
	prod_dnorm_quant<-dnorm_quant*quant_norm
	dnorm_prod_diff<-c(prod_dnorm_quant,0)-c(0,prod_dnorm_quant)
	n_interv_adj<-n_interv-1
	negscore_mu<-(-sum(quant_norm)+sum(dnorm_diff/pnorm_diff*n_interv_adj))/sigma
	negscore_sigma<-(-sum(quant_norm^2)+d+sum(dnorm_prod_diff/pnorm_diff*n_interv_adj))/sigma
	neglikscore_val<-c(negscore_mu,negscore_sigma)
	return(neglikscore_val)
}
prior<-function(theta){ # Uniform prior density 
	mu<-theta[1]; sigma<-abs(theta[2])
	if(mu>mu_max||mu<mu_min) return(0)
	if(sigma>sigma_max||sigma<sigma_min) return(0)
	return(1)
}

	<-function(theta,sumstat){ # Output the logarithm of the unnormalized posterior density
	loglik_val<--negloglik(theta,sumstat)
	prior_val<-prior(theta)
	if(prior_val==1) logprior_val<-0
	if(prior_val==0) logprior_val<--Inf
	return(loglik_val+logprior_val)
}

dist_mat_vec<-function(A,v){ # Output the Euclidean distane between each column of matrix A and vector v
	distance<-c(sqrt(colSums((A-v)^2)))
	return(distance)
}

accept_kernel<-function(v){ # Input the distance between summaries scaled by the bandwidth; Output decisions of acceptance/rejection
	N0<-length(v)
	acc_prob<-exp(-v^2*pi) # Use N(0,1/(2pi)) kernel density
	acc_unif<-runif(N0) # Generate uniform r.v. for accept/reject
	acc_ind<-acc_unif<acc_prob
	return(acc_ind)
}

diff_acc_rate<-function(tol,v,want_acc_rate){# Used in optimization to find the tolerance giving a certain acceptance rate; Input the tolerance, distance between summaries and desired accepted rate; Ouput the difference between the expected and desired accepted rates
	scaled_v<-v/tol
	acc_prob<-exp(-scaled_v^2*pi) # Use N(0,1/(2pi)) kernel density
	exp_acc_rate<-sum(acc_prob)/length(v)
	return(abs(exp_acc_rate-want_acc_rate))
}

sign_sq<-function(x) return(sign(x)*x^2)

CalSummary<-function(data,sum_pars){ # Output the quantiles of the dataset and their exponentials
	probs<-sum_pars
	quant<-quantile(data,probs=probs)
	sumstat<-sign_sq(quant)
	# sumstat<-exp(quant/2)
	# sumstat<-quant
	return(sumstat)
}

summary_simu<-function(MC_theta,n,sum_pars){ # Output the simulated summary statistics from the pre-generated pool of normal quantiles
	probs_i<-match(length(sum_pars),probs_size_all)
	n_i<-match(n,n_all); d<-length(sum_pars)
	MC_mu<-MC_theta[,1]; MC_sigma<-MC_theta[,2]; N0<-length(MC_mu)
	base_quant<-normquantile_pool[[n_i]][1:N0,probs_indices_pool[[probs_i]]] # N0*d
	quant<-base_quant*MC_sigma+MC_mu # N0*d
	sumstat<-sign_sq(quant)
	#sumstat<-exp(quant/2)
	# sumstat<-quant
	return(sumstat)
}

summary_dimreduce_matrix<-function(sum_pars){ # Output the dimension reduction matrix based on the analytical for of the asymptotic distribution of the summaries
	probs<-sum_pars; probs_size<-length(probs)
	stanardnorm_quant<-qnorm(probs)
	dens_quant<-dnorm(stanardnorm_quant)
	asy_var_part<-matrix(0,probs_size,probs_size)
	for(i in 1:probs_size){for(j in 1:probs_size) asy_var_part[i,j]<-probs[min(i,j)]*(1-probs[max(i,j)])}
	asy_var<-diag(1/dens_quant)%*%asy_var_part%*%diag(1/dens_quant)
	Ds<-cbind(rep(1,probs_size),stanardnorm_quant)
	dimreduc_matrix<-t(Ds)%*%solve(asy_var)/sigma_targ^2
	#dimreduc_matrix_exp<-4*dimreduc_matrix%*%diag(exp(-sigma_targ*stanardnorm_quant/2-mu_targ/2))
	dimreduc_matrix_sgnsq<-dimreduc_matrix%*%diag(abs(sigma_targ*stanardnorm_quant+mu_targ)^(-1))
	return(dimreduc_matrix_sgnsq)
}

IS_ABC<-function(N,prop_mean,prop_sd,n,sum_pars,sumstat_obs,reduce_dim=FALSE,acc_rate=0){
# This function is used to pre-calculate the tolerance value for certain proposal distribution and acceptance rate
# Output the tolerance value
# N is the MC size; Uniform proposal distribution is used so that the importance weights are 1 or 0; prop_mean and prop_sd are the mean and sd of the uniform proposal; n is the data size; sum_pars is the probabilities for the quantile summaries; sumstat_obs is the observed summary; reduce_dim indicates whether or not reducing the dimension
		if(!reduce_dim){
			MC_theta<-t(rbind(runif(N,-1,1),runif(N,-1,1))*sqrt(3)*prop_sd+prop_mean) # N*p; uniform proposal distribution
			MC_sumstat<-summary_simu(MC_theta,n,sum_pars) # N*d
			MC_dist<-dist_mat_vec(t(MC_sumstat),sumstat_obs) # N vector
			quantile_tol_results<-optimise(diff_acc_rate,interval=c(min(MC_dist),max(MC_dist)),v=MC_dist,want_acc_rate=acc_rate)
			quantile_tol<-quantile_tol_results$minimum
		}
		if(reduce_dim){
			probs_i<-match(length(sum_pars),probs_size_all)
			sumstat_obs_lincomb<-c(linear_comb_matrix[[probs_i]]%*%sumstat_obs) # 2*1
			MC_theta<-t(rbind(runif(N,-1,1),runif(N,-1,1))*sqrt(3)*prop_sd+prop_mean) # N*p
			MC_sumstat<-summary_simu(MC_theta,n,sum_pars) # N*d
			MC_sumstat_lincomb<-linear_comb_matrix[[probs_i]]%*%t(MC_sumstat) # 2*N
			MC_dist<-dist_mat_vec(MC_sumstat_lincomb,sumstat_obs_lincomb) # N vector
			quantile_tol_results<-optimise(diff_acc_rate,interval=c(min(MC_dist),max(MC_dist)),v=MC_dist,want_acc_rate=acc_rate)
			quantile_tol<-quantile_tol_results$minimum
		}
		return(list(quant_tol=quantile_tol))
}

IS_ABC_fixedMC_fixedtol<-function(N,prop_mean,prop_sd,n,sum_pars,sumstat_obs,reduce_dim=FALSE,tol=1,output_tol_i=0){
# This function is the ABC algorithm using a fixed MC size and a series of pre-calculated tolerance values
# Output the estimated theta and estimated IS variance
# Input parameters are the same as the above; output_tol_i indicates whether or not output the simulated and accepted parameter values and summaries
	tol_length<-length(tol)
	theta_est<-matrix(0,tol_length,2); thetasqr_est<-theta_est
	theta_est_regadjused<-theta_est; thetasqr_est_regadjused<-theta_est
	postvar_est<-theta_est; postvar_est_regadjused<-theta_est
	output<-FALSE
	if(!reduce_dim){
		MC_theta<-t(rbind(runif(N,-1,1),runif(N,-1,1))*sqrt(3)*prop_sd+prop_mean) # N*p; uniform proposal distribution
		MC_sumstat<-summary_simu(MC_theta,n,sum_pars) # N*d
		MC_dist<-dist_mat_vec(t(MC_sumstat),sumstat_obs) # N vector
		for(tol_i in 1:tol_length){
			tolerance<-tol[tol_i]
			MC_dist_scaled<-MC_dist/tolerance
			acc_ind<-accept_kernel(MC_dist_scaled); N_acc<-sum(acc_ind)
			theta_acc<-MC_theta[acc_ind,]
			sumstat_acc<-MC_sumstat[acc_ind,]

			# Basic variance estimator
			theta_est[tol_i,]<-colMeans(theta_acc)
			thetasqr_est[tol_i,]<-colMeans(theta_acc^2)
			postvar_est[tol_i,]<-thetasqr_est[tol_i,]-theta_est[tol_i,]^2

			# Regression-adjusted variance estimator
			regresults<-lm(theta_acc~sumstat_acc)
			regcoef<-regresults$coef
#browser()
			theta_acc_regadjused<-theta_acc-sumstat_acc%*%regcoef[-1,]
			theta_est_regadjused[tol_i,]<-colMeans(theta_acc_regadjused)
			thetasqr_est_regadjused[tol_i,]<-colMeans(theta_acc_regadjused^2)
			postvar_est_regadjused[tol_i,]<-thetasqr_est_regadjused[tol_i,]-theta_est_regadjused[tol_i,]^2

			if(tol_i==output_tol_i){
				output<-TRUE 
				MC_theta_output<-MC_theta 
				theta_acc_output<-theta_acc
				output_results<-list(MC=MC_theta_output,acc=theta_acc_output)
			}
		}
	}
	if(reduce_dim){
		probs_i<-match(length(sum_pars),probs_size_all)
		sumstat_obs_lincomb<-c(linear_comb_matrix[[probs_i]]%*%sumstat_obs) # 2*1
		MC_theta<-t(rbind(runif(N,-1,1),runif(N,-1,1))*sqrt(3)*prop_sd+prop_mean) # N*p
		MC_sumstat<-summary_simu(MC_theta,n,sum_pars) # N*d
		MC_sumstat_lincomb<-linear_comb_matrix[[probs_i]]%*%t(MC_sumstat) # 2*N
		MC_dist<-dist_mat_vec(MC_sumstat_lincomb,sumstat_obs_lincomb) # N vector
		for(tol_i in 1:tol_length){
			tolerance<-tol[tol_i]
			MC_dist_scaled<-MC_dist/tolerance
			acc_ind<-accept_kernel(MC_dist_scaled); N_acc<-sum(acc_ind)
			theta_acc<-MC_theta[acc_ind,]
			sumstat_acc<-MC_sumstat[acc_ind,]

			# Basic variance estimator
			theta_est[tol_i,]<-colMeans(theta_acc)
			thetasqr_est[tol_i,]<-colMeans(theta_acc^2)
			postvar_est[tol_i,]<-thetasqr_est[tol_i,]-theta_est[tol_i,]^2

			# Regression-adjusted variance estimator
			regresults<-lm(theta_acc~sumstat_acc)
			regcoef<-regresults$coef
			theta_acc_regadjused<-theta_acc-sumstat_acc%*%regcoef
			theta_est_regadjused[tol_i,]<-colMeans(theta_acc_regadjused)
			thetasqr_est_regadjused[tol_i,]<-colMeans(theta_acc_regadjused^2)
			postvar_est_regadjused[tol_i,]<-thetasqr_est_regadjused[tol_i,]-theta_est_regadjused[tol_i,]^2

			if(tol_i==output_tol_i){
				output<-TRUE 
				MC_theta_output<-MC_theta 
				theta_acc_output<-theta_acc
				output_results<-list(MC=MC_theta_output,acc=theta_acc_output)
			}
		}
	}
	if(output) return(list(postvar=postvar_est,postvar_reg=postvar_est_regadjused,output_results=output_results))
	return(list(postvar=postvar_est,postvar_reg=postvar_est_regadjused))
}

#########################################################################
# Finished experiments and datasets
#########################################################################
# The experiments take a long time. Here are some finished results in the workspace files.


# This contains simulated normal quantiles for n=10^(2:5) and probs_size_all=c(2,4,9,19); 
# normquantile_pool is a 1e6*sum(probs_size_all) matrix containing the simulated quantiles
# probs_indices_pool contains the indicies of each set of probabilities for tested quantiles
load("E:\\normquantile_pool2.RData")

# Experiment for n=1e2, replicated  times
load("E:\\N1e6_experiment.RData")
# Experiment for n=1e3, replicated  times
# Experiment for n=1e4, replicated  times
# Experiment for n=1e5, replicated  times

save.image('E:\\n1e2_experiment.RData')
save.image('E:\\n1e3_experiment.RData')
save.image('E:\\n1e4_experiment.RData')
save.image('E:\\n1e5_experiment.RData')


load('E:\\n1e2_experiment.RData')
load('E:\\n1e3_experiment.RData')
load('E:\\n1e4_experiment.RData')
load('E:\\n1e5_experiment.RData')



#######################################################################
# Parameter settings
#######################################################################
random_seed<-100
mu_targ<-1; sigma_targ<-sqrt(2); theta_targ<-c(mu_targ,sigma_targ)
mu_min<--10; mu_max<-10; sigma_min<--10; sigma_max<-10
n_all<-10^(2:5); n_max<-max(n_all); n_length<-length(n_all)
probs_size_all<-c(2,4,9,19); probs_length<-length(probs_size_all)

probs_all<-list(0) # Calculate the tested quantile probabilities
for(probs_i in 1:probs_length) {
	probs_size<-probs_size_all[probs_i]; 
	probs_all[[probs_i]]<-seq(0,1,length=probs_size_all[probs_i]+2)[c(-1,-(probs_size+2))]
}


#########################################################################################################################
#########################################################################################################################
# Parameter estimations using MLE and MCMC 
#########################################################################################################################
#########################################################################################################################
replic<-100
theta_MLES_mu_all<-array(0,c(replic,probs_length,n_length)); theta_MLES_sigma_all<-array(0,c(replic,probs_length,n_length))
theta_MLE_mu_all<-array(0,c(replic,n_length)); theta_MLE_sigma_all<-array(0,c(replic,n_length))
theta_MCMC_mu_all<-array(0,c(replic,probs_length,n_length)); theta_MCMC_sigma_all<-array(0,c(replic,probs_length,n_length))
postvar_MCMC_mu_all<-array(0,c(replic,probs_length,n_length)); postvar_MCMC_sigma_all<-array(0,c(replic,probs_length,n_length))

set.seed(100); plot(0,0,xlim=c(1,replic),ylim=c(1,replic))
for(rep_i in 1:replic){ # Loop of experiment replication with different observations
	points(rep_i,rep_i,pch=20)
	for(n_i in 1:n_length){ # Loop of different n
		n<-n_all[n_i]; observ<-rnorm(n,mu_targ,sigma_targ)
		theta_MLE_mu_all[rep_i,n_i]<-mean(observ); theta_MLE_sigma_all[rep_i,n_i]<-sd(observ)
		for(probs_i in 1:probs_length){ # Loop of different choices of summary
			probs<-probs_all[[probs_i]]; d<-length(probs)
			summary_observ<-CalSummary(observ,probs)
			n_interv<-c(floor(n*probs),n)-c(0,floor(n*probs))
			
			# MLES
			theta_MLES_results<-optim(theta_targ,negloglik,gr=neglikscore,method='BFGS',control=list(trace=0,fnscale=n),sumstat=summary_observ)
			theta_MLES<-theta_MLES_results$par
			theta_MLES_mu_all[rep_i,probs_i,n_i]<-theta_MLES[1]
			theta_MLES_sigma_all[rep_i,probs_i,n_i]<-theta_MLES[2]

			# Posterior mean by MCMC
			MCMC_results<-metrop(logposterior,initial=theta_targ,scale=diag(0.03,2),nbatch=1e4,  sumstat=summary_observ)
			theta_MCMC<-colMeans(MCMC_results$batch)
			thetasqr_MCMC<-colMeans(MCMC_results$batch^2)
			postvar_MCMC<-thetasqr_MCMC-theta_MCMC^2			
			postvar_MCMC_mu_all[rep_i,probs_i,n_i]<-postvar_MCMC[1]
			postvar_MCMC_sigma_all[rep_i,probs_i,n_i]<-postvar_MCMC[2]
		}
	}
}
dev.off()

# Calculating MSE of MLES and MLE
MSE_MLES<-array(0,c(n_length,probs_length,2)); MSE_MLE<-matrix(0,n_length,2)
MSE_MCMC<-array(0,c(n_length,probs_length,2))
for(n_i in 1:n_length){
	n<-n_all[n_i]
	for(probs_i in 1:probs_length){
		MSE_MLES[n_i,probs_i,1]<-mean((theta_MLES_mu_all[,probs_i,n_i]-mu_targ)^2)*n
		MSE_MLES[n_i,probs_i,2]<-mean((theta_MLES_sigma_all[,probs_i,n_i]-sigma_targ)^2)*n
		MSE_MCMC[n_i,probs_i,1]<-mean((theta_MCMC_mu_all[,probs_i,n_i]-mu_targ)^2)*n
		MSE_MCMC[n_i,probs_i,2]<-mean((theta_MCMC_sigma_all[,probs_i,n_i]-sigma_targ)^2)*n	
	}
	MSE_MLE[n_i,1]<-mean((theta_MLE_mu_all[,n_i]-mu_targ)^2)*n
	MSE_MLE[n_i,2]<-mean((theta_MLE_sigma_all[,n_i]-sigma_targ)^2)*n			
}

################################################################################################
################################################################################################
# Parameter estimations using ABC
################################################################################################
################################################################################################

#################################################################################################
# Generate observations and their summaries
#################################################################################################
#MLE and MCMC (smaller replication for reasonable computing time)
replic<-100; N<-1e5
# Data structure initialization
theta_MLES_mu_all<-array(0,c(replic,probs_length,n_length)); theta_MLES_sigma_all<-array(0,c(replic,probs_length,n_length))
theta_MCMC_mu_all<-array(0,c(replic,probs_length,n_length)); theta_MCMC_sigma_all<-array(0,c(replic,probs_length,n_length))
theta_MLE_mu_all<-array(0,c(replic,n_length)); theta_MLE_sigma_all<-array(0,c(replic,n_length))

summary_observ_all<-list(0) # For saving the observed summaries
for(probs_i in 1:probs_length) {
	probs_size<-probs_size_all[probs_i]; 
	summary_observ_all[[probs_i]]<-array(0,c(replic,n_length,probs_size))
}

# This part is the same as the part of MLES estimation
set.seed(100); plot(0,0,xlim=c(1,replic),ylim=c(1,replic))
for(rep_i in 1:replic){ # Loop of experiment replication with different observations
	observ_all<-rnorm(n_max,mu_targ,sigma_targ)
	points(rep_i,rep_i,pch=20)
	for(n_i in 1:n_length){
		n<-n_all[n_i]; observ<-observ_all[1:n]		
		theta_MLE_mu_all[rep_i,n_i]<-mean(observ); theta_MLE_sigma_all[rep_i,n_i]<-sd(observ)

		for(probs_i in 1:probs_length){ # Loop of different choices of summary
			probs<-probs_all[[probs_i]]; d<-length(probs)
			summary_observ<-CalSummary(observ,probs)
			summary_observ_all[[probs_i]][rep_i,n_i,]<-summary_observ
			n_interv<-c(floor(n*probs),n)-c(0,floor(n*probs))
			
			# MLES, in order to obtain the distribution of the MLES 
			theta_MLES_results<-optim(theta_targ,negloglik,gr=neglikscore,method='BFGS',control=list(trace=0,fnscale=n),sumstat=summary_observ)
			theta_MLES<-theta_MLES_results$par
			theta_MLES_mu_all[rep_i,probs_i,n_i]<-theta_MLES[1]
			theta_MLES_sigma_all[rep_i,probs_i,n_i]<-theta_MLES[2]
		}
	}
}
dev.off()


#######################################################
# Calculate the dimension reduction matrices
#######################################################
linear_comb_matrix<-list(0) 
for(probs_i in 1:probs_length){
	probs<-probs_all[[probs_i]]
	linear_comb_matrix[[probs_i]]<-summary_dimreduce_matrix(probs)
}

#####################################################################
# Proposal distribution setting; Uniform is used as the proposals
#####################################################################
# Use the sd of the MLES estimator as the base sd of the proposal distributions
set.seed(random_seed)
prop_sd_all<-array(0,c(probs_length,n_length,2))
for(n_i in 1:n_length){
	for(probs_i in 1:probs_length){
		prop_sd_all[probs_i,n_i,]<-c(sd(theta_MLES_mu_all[,probs_i,n_i]),sd(theta_MLES_sigma_all[,probs_i,n_i]))
	}
}

# Construct the proposal means by adding a constant, calculated by multiplying mean_fac with the base sd above, to the MLES estimator, and the proposal sd by scaling some power of the base sd with sd_fac;
no_fac<-2; no_power<-3; no_prop<-no_fac*no_power
mean_fac_all<-c(.25,.5,.5,1,.5,1)
sd_fac_all<-6*mean_fac_all
power_all<-rep(c(1,2/3,1/2),c(2,2,2)); 
propal_tune_all<-cbind(mean_fac_all,sd_fac_all,power_all)
dimnames(propal_tune_all)=list(NULL,c("mean_fac","sd_fac","power"))

prop_pars_adj_all<-list(0)
for(n_i in 1:n_length){
	prop_pars_adj_all[[n_i]]<-list(0)
	for(probs_i in 1:probs_length){
		prop_pars_adj_all[[n_i]][[probs_i]]<-list(mean=matrix(0,no_prop,2),sd=matrix(0,no_prop,2))
		mean_dev<-prop_sd_all[probs_i,n_i,]
		prop_mean_adj<-cbind(propal_tune_all[,1]*mean_dev[1]^propal_tune_all[,3],propal_tune_all[,1]*mean_dev[2]^propal_tune_all[,3])
		prop_sd_adj<-cbind(propal_tune_all[,2]*mean_dev[1]^propal_tune_all[,3],propal_tune_all[,2]*mean_dev[2]^propal_tune_all[,3])
		prop_pars_adj_all[[n_i]][[probs_i]]$mean<-prop_mean_adj
		prop_pars_adj_all[[n_i]][[probs_i]]$sd<-prop_sd_adj
	}
}

#######################################################
# Data structure initialization
#######################################################
# Since it is difficult to test a wide enough range of tolerance with one proposal distribution, 6
# proposal distributions with different means and variances are tested for the same series of
# acceptance rates and the estimated theta are recorded. Then for each choice of summary, calculate
# MSE for each setting of (prop_mean, prop_sd, acceptance rate) and combine together in the order of
# tolerance values. 

series_length<-100; acc_rate_series<-seq(0,.1,length=series_length+1)[-1]
ttlseries_length<-series_length*no_prop
probs_lincomb_length<-probs_length*2 # (probs_length+1):probs_lincomb_length are for the results of dimension reduced summaries
probs_lin_names<-c(probs_size_all[1:probs_length],paste(probs_size_all[1:probs_length],'_lin',sep=''))
tol_all<-array(0,c(no_prop,series_length,probs_lincomb_length,n_length))
dimnames(tol_all)<-list(NULL,1:series_length,probs_lin_names,log10(n_all)[1:n_length])
postvar_ABC_all<-rep(list(array(0,c(no_prop,series_length,replic,probs_lincomb_length,n_length))),2)
dimnames(postvar_ABC_all[[1]])<-list(NULL,1:series_length,NULL,probs_lin_names,log10(n_all)[1:n_length])
dimnames(postvar_ABC_all[[2]])<-dimnames(postvar_ABC_all[[1]])
ISvar_ABC_all<-postvar_ABC_all

#######################################################
# Pre-calculate the series of tolerance levels
#######################################################

n_simu<-4 # Tested data size
set.seed(random_seed); rep_i<-1
x11(); plot(0,0,xlim=c(1,n_length*probs_length),ylim=c(1,series_length))
#quartz(); plot(0,0,xlim=c(1,n_length*probs_length),ylim=c(1,series_length))
for(n_i in n_simu){ # Loop of data size
	n<-n_all[n_i]
	for(probs_i in 1:probs_length){ # Loop of choices of summary
		probs<-probs_all[[probs_i]]
		summary_observ<-summary_observ_all[[probs_i]][rep_i,n_i,]
		theta_MLES<-c(theta_MLES_mu_all[rep_i,probs_i,n_i],theta_MLES_sigma_all[rep_i,probs_i,n_i])
		mean_dev<-prop_sd_all[probs_i,n_i,]; 

		for(acc_i in 1:series_length){ # Loop of acceptance rate
			points((n_i-1)*probs_length+probs_i,acc_i,pch=20)
			acc_rate<-acc_rate_series[acc_i]

			for(prop_i in 1:no_prop){ # Loop of proposal distribution
				prop_mean<-theta_MLES-prop_pars_adj_all[[n_i]][[probs_i]]$mean[prop_i,]
				prop_sd<-prop_pars_adj_all[[n_i]][[probs_i]]$sd[prop_i,]
				ABC_results<-IS_ABC(N,prop_mean,prop_sd,n,probs,summary_observ,acc_rate=acc_rate)
				ABC_lin_results<-IS_ABC(N,prop_mean,prop_sd,n,probs,summary_observ,reduce_dim=T,acc_rate=acc_rate)
				tol<-ABC_results$quant_tol; tol_lin<-ABC_lin_results$quant_tol
				tol_all[prop_i,acc_i,probs_i,n_i]<-tol
				tol_all[prop_i,acc_i,probs_i+probs_length,n_i]<-tol_lin
			}
			
		}
	}
}
dev.off()

#######################################################
# ABC estimators for the series of tolerance levels
#######################################################
set.seed(random_seed)
#x11(); plot(0,0,xlim=c(1,replic),ylim=c(1,replic))
quartz(); plot(0,0,xlim=c(1,replic),ylim=c(1,replic))
for(rep_i in 1:replic){ # Loop of experiment replication with different observations
	points(rep_i,rep_i,pch=20)
	for(n_i in n_simu){ # Loop of data size
		n<-n_all[n_i]
		for(probs_i in 1:probs_length){ # Loop of different choices of summary
			probs<-probs_all[[probs_i]]
			summary_observ<-summary_observ_all[[probs_i]][rep_i,n_i,]
			theta_MLES<-c(theta_MLES_mu_all[rep_i,probs_i,n_i],theta_MLES_sigma_all[rep_i,probs_i,n_i])
			mean_dev<-prop_sd_all[probs_i,n_i,]; 

			# ABC posterior mean
			for(prop_i in 1:no_prop){
				prop_mean<-theta_MLES-prop_pars_adj_all[[n_i]][[probs_i]]$mean[prop_i,]
				prop_sd<-prop_pars_adj_all[[n_i]][[probs_i]]$sd[prop_i,]
				ABC_results<-IS_ABC_fixedMC_fixedtol(N,prop_mean,prop_sd,n,probs,summary_observ,tol=tol_all[prop_i,,probs_i,n_i])

				postvar_ABC_all[[1]][prop_i,,rep_i,probs_i,n_i]<-ABC_results$postvar[,1]
				postvar_ABC_all[[2]][prop_i,,rep_i,probs_i,n_i]<-ABC_results$postvar[,2]
				postvar_ABC_all[[1]][prop_i,,rep_i,probs_i+probs_length,n_i]<-ABC_results$postvar_reg[,1]
				postvar_ABC_all[[2]][prop_i,,rep_i,probs_i+probs_length,n_i]<-ABC_results$postvar_reg[,2]
			}
		}
	}
}
dev.off()

#######################################################
# Plots of tolerances vs MSE
# The performance of sigma estimator with reduced-dimension
# statistic is much better than that with original. But 
# the bias of the mu estimator is larger and larger than that 
# of the original as the tolerance increases, which may be due to 
# that the second derivatives are increased by the transformation
# matrix.
#######################################################
plot_tol_length<-50; plot_prop_length<-4
plot_length<-ttlseries_length/series_length*plot_tol_length/no_prop*plot_prop_length # length of the tolerances selected to be drawn in the plot
MSE_all<-array(0,c(plot_length,n_length,probs_lincomb_length,2),dimnames=list(NULL,log10(n_all)[1:n_length],probs_lin_names,c('mu','sigma')))
ISvar_all<-MSE_all; var_all<-MSE_all; tol_plot_all<-MSE_all
MSE_plot<-MSE_all; ISvar_plot<-MSE_all; var_plot<-MSE_all; tol_plot<-MSE_all

for(est_par_i in 1:2){
	for(n_i in n_simu){
		n<-n_all[n_i]
		for(probs_i in 1:probs_length){
			plot_i<-0
			for(ttlseries_i in 1:ttlseries_length){
				series_i<-floor((ttlseries_i-1)/no_prop)+1 # indices of the acceptance rates
				prop_i<-ttlseries_i-(series_i-1)*no_prop # indices of the proposal distributions
				if(prop_i<=plot_prop_length&series_i<=plot_tol_length){ # This controls the output of which proposal distributions and which acceptance rates should be drawn
					plot_i<-plot_i+1
						MSE_all[plot_i,n_i,probs_i,est_par_i]<-mean((theta_ABC_all[[est_par_i]][prop_i,series_i,,probs_i,n_i]-theta_targ[est_par_i])^2)*n
						var_all[plot_i,n_i,probs_i,est_par_i]<-var(theta_ABC_all[[est_par_i]][prop_i,series_i,,probs_i,n_i])*n	
						ISvar_all[plot_i,n_i,probs_i,est_par_i]<-mean(ISvar_ABC_all[[est_par_i]][prop_i,series_i,,probs_i,n_i])*n
						MSE_all[plot_i,n_i,probs_i+probs_length,est_par_i]<-mean((theta_ABC_all[[est_par_i]][prop_i,series_i,,probs_i+probs_length,n_i]-theta_targ[est_par_i])^2)*n
						var_all[plot_i,n_i,probs_i+probs_length,est_par_i]<-var(theta_ABC_all[[est_par_i]][prop_i,series_i,,probs_i+probs_length,n_i])*n	
						ISvar_all[plot_i,n_i,probs_i+probs_length,est_par_i]<-mean(ISvar_ABC_all[[est_par_i]][prop_i,series_i,,probs_i+probs_length,n_i])*n
						tol_plot_all[plot_i,n_i,probs_i,est_par_i]<-tol_all[prop_i,series_i,probs_i,n_i]	
				}
			}
			order_tol<-order(tol_plot_all[,n_i,probs_i,est_par_i]) # sort the tolerance values
			MSE_plot[,n_i,probs_i,est_par_i]<-MSE_all[order_tol,n_i,probs_i,est_par_i]
			ISvar_plot[,n_i,probs_i,est_par_i]<-ISvar_all[order_tol,n_i,probs_i,est_par_i]
			var_plot[,n_i,probs_i,est_par_i]<-var_all[order_tol,n_i,probs_i,est_par_i]
			MSE_plot[,n_i,probs_i+probs_length,est_par_i]<-MSE_all[order_tol,n_i,probs_i+probs_length,est_par_i]
			ISvar_plot[,n_i,probs_i+probs_length,est_par_i]<-ISvar_all[order_tol,n_i,probs_i+probs_length,est_par_i]
			var_plot[,n_i,probs_i+probs_length,est_par_i]<-var_all[order_tol,n_i,probs_i+probs_length,est_par_i]
			tol_plot[,n_i,probs_i,est_par_i]<-tol_plot_all[order_tol,n_i,probs_i,est_par_i]
		}
	}
}
MSE_hABC_all<-MSE_plot-ISvar_plot

x11(); par(mfcol=c(4,2))
#quartz(); par(mfrow=c(4,2))
for(est_par_i in 1:2){
	for(n_i in n_simu){
		n<-n_all[n_i]
		ymax<-max(c(MSE_plot[-1,n_i,,est_par_i],MSE_MLES[n_i,,est_par_i],MSE_MLE[n_i,est_par_i]))
		ymin<-min(c(MSE_plot[-1,n_i,,est_par_i],MSE_MLES[n_i,,est_par_i],MSE_MLE[n_i,est_par_i]))
		for(probs_i in 1:probs_length){
			xmax<-max(c(tol_plot[,n_i,probs_i,est_par_i],tol_plot[,n_i,probs_i+probs_length,est_par_i]))
			xmin<-min(c(tol_plot[,n_i,probs_i,est_par_i],tol_plot[,n_i,probs_i+probs_length,est_par_i]))
	
			# Some jigsaw shape may be observed in some parts of the curves in the plot. Because one issue of combining the output from different proposal distributions is that, for the same values of the tolerance level, different proposal distributions could give different MSE. The reason of this issue is the importance sampling bias for larger acceptance rate. Therefore if restrict the acceptance rate under 5% or use larger MC size, this issue will be mitigated.
			# plot(tol_plot[,n_i,probs_i],MSE_plot[,n_i,probs_i],type='l',main=paste('n=',n_all[n_i],', summary=',probs_size_all[probs_i],sep=''),ylim=c(0,ymax),ylab='n*MSE',xlab=bquote(epsilon))
			# lines(tol_plot[,n_i,probs_i],MSE_plot[,n_i,probs_i+probs_length],lty=2)
			# abline(h=MSE_MLES[n_i,probs_i,est_par_i],col=2)
			# abline(h=MSE_MLE[n_i,est_par_i],col=3)
			# This is to make a smoother curve, averaging the MSE from different proposal distributions. 
			plot(ksmooth(tol_plot[,n_i,probs_i,est_par_i],MSE_plot[,n_i,probs_i,est_par_i],kernel='box',bandwidth=0.01),type='l',main=paste('n=',n_all[n_i],', summary=',probs_size_all[probs_i],sep=''),ylim=c(ymin-0.05,ymax),ylab='MSE',xlab=bquote(epsilon))
			lines(ksmooth(tol_plot[,n_i,probs_i,est_par_i],MSE_plot[,n_i,probs_i+probs_length,est_par_i],kernel='box',bandwidth=0.01),lty=2)
			#lines(ksmooth(tol_plot[,n_i,probs_i],var_plot[,n_i,probs_i],kernel='box',bandwidth=0.1),col=2)
			#lines(ksmooth(tol_plot[,n_i,probs_i],var_plot[,n_i,probs_i+probs_length],kernel='box',bandwidth=0.1),col=2,lty=2)
			abline(h=MSE_MLES[n_i,probs_i,est_par_i],col=2)
			abline(h=MSE_MLE[n_i,est_par_i],col=3)
		}
	}
}

quartz(); par(mfrow=c(2,2))
for(n_i in n_simu){
	n<-n_all[n_i]
	for(probs_i in 1:probs_length){
		ymax<-max(c(MSE_plot[-1,n_i,probs_i],MSE_plot[-1,n_i,probs_i+probs_length],MSE_MLES[n_i,probs_i],MSE_MLE[n_i]))
		ymin<-min(c(MSE_plot[-1,n_i,probs_i],MSE_plot[-1,n_i,probs_i+probs_length],MSE_MLES[n_i,probs_i],MSE_MLE[n_i]))
		xmax<-max(c(tol_plot[,n_i,probs_i],tol_plot[,n_i,probs_i+probs_length]))
		xmin<-min(c(tol_plot[,n_i,probs_i],tol_plot[,n_i,probs_i+probs_length]))
		if(probs_i>1) title<-paste(probs_size_all[probs_i],' summaries',sep='')
		if(probs_i==1) title<-paste(1,' summary',sep='')

		plot(ksmooth(tol_plot[,n_i,probs_i],MSE_plot[,n_i,probs_i],kernel='box',bandwidth=0.1),type='l',main=title,ylim=c(0,15),ylab='n*MSE',xlab=bquote(epsilon))
		lines(ksmooth(tol_plot[,n_i,probs_i],MSE_plot[,n_i,probs_i+probs_length],kernel='box',bandwidth=0.1),lty=2)
		abline(h=MSE_MLES[n_i,probs_i],col=2)
		abline(h=MSE_MLE[n_i],col=3)
		if(probs_i==1) legend("topleft",c("Original","Dimension Reduced","MLES","MLE"),lty=c(1,2,1,1),col=c(1,1,2,3))
	}
}

quartz(); par(mfrow=c(2,2))
for(n_i in n_simu){
	n<-n_all[n_i]
	for(probs_i in 1:probs_length){
		ymax<-max(c(var_plot[-1,n_i,probs_i],var_plot[-1,n_i,probs_i+probs_length],MSE_MLES[n_i,probs_i],MSE_MLE[n_i]))
		ymin<-min(c(var_plot[-1,n_i,probs_i],var_plot[-1,n_i,probs_i+probs_length],MSE_MLES[n_i,probs_i],MSE_MLE[n_i]))
		xmax<-max(c(tol_plot[,n_i,probs_i],tol_plot[,n_i,probs_i+probs_length]))
		xmin<-min(c(tol_plot[,n_i,probs_i],tol_plot[,n_i,probs_i+probs_length]))

		plot(ksmooth(tol_plot[,n_i,probs_i],var_plot[,n_i,probs_i],kernel='box',bandwidth=0.1),type='l',main=paste(probs_size_all[probs_i],' summaries',sep=''),ylim=c(0,7),ylab='n*MSE',xlab=bquote(epsilon))
		lines(ksmooth(tol_plot[,n_i,probs_i],var_plot[,n_i,probs_i+probs_length],kernel='box',bandwidth=0.1),lty=2)
		abline(h=MSE_MLES[n_i,probs_i],col=2)
		abline(h=MSE_MLE[n_i],col=3)
		if(probs_i==1) legend("topleft",c("Original","Dimension Reduced","MLES","MLE"),lty=c(1,2,1,1),col=c(1,1,2,3))
	}
}


#x11(); par(mfrow=c(2,2))
quartz(); par(mfrow=c(2,2))
for(n_i in n_simu){
	n<-n_all[n_i]
	for(probs_i in 1:probs_length){
		ymax<-max(c(ISvar_plot[,n_i,probs_i],ISvar_plot[,n_i,probs_i+probs_length]))
		ymin<-min(c(ISvar_plot[,n_i,probs_i],ISvar_plot[,n_i,probs_i+probs_length]))
		xmax<-max(c(tol_plot[,n_i,probs_i],tol_plot[,n_i,probs_i+probs_length]))
		xmin<-min(c(tol_plot[,n_i,probs_i],tol_plot[,n_i,probs_i+probs_length]))
		plot(tol_plot[,n_i,probs_i],ISvar_plot[,n_i,probs_i],type='l',main=paste('n=',n_all[n_i],', summary=',probs_size_all[probs_i],sep=''),ylim=c(ymin,ymax))
		lines(tol_plot[,n_i,probs_i],ISvar_plot[,n_i,probs_i+probs_length],lty=2)
		#abline(h=MSE_MLES[n_i,probs_i],col=2)
		#abline(h=MSE_MLE[n_i],col=3)
	}
}









