These are the codes for the paper “Convergence of Regression Adjusted Approximate Bayesian
Computation”. Biometrika, 105(2):301-318.

## Dependencies

The scripts rely on the following R packages:

  * `glmnet`
  * `matrixStats`
  * `ks`
  * `pbapply`
  * `snowfall`
  * `parallel`
  * `MASS`
  * `mvtnorm`
  * `tseries`
  * `numDeriv`
  * `mcmc`
-----

## File Summaries

The repository uses the 'Utility functions' repository. These scripts are numerical studies that apply the utility functions to specific statistical models, as seen in the paper's numerical examples.

  * **`Regression Gaussian example.R`**: An analysis script for a univariate Gaussian (normal) model. It serves as a benchmark to compare the performance (Mean Squared Error) of standard MLE, MCMC, standard ABC, and the paper's primary topic: **Regression-Adjusted ABC**.

  * **`SVAR_1.R`**: A detailed analysis of a Stochastic Volatility with AR(1) dynamics (SVAR(1)) model. This script compares two advanced estimation methods:

    1.  **Efficient Method of Moments (EMM)**: Uses the score of an auxiliary GARCH(1,1) model as summary statistics.
    2.  **Indirect Inference (II)**: Uses simpler moment-based summaries (mean, variance, autocovariance).

    For both methods, the script demonstrates a "Reverse Sampler" and uses the estimates to build an efficient proposal distribution for **Regression Importance-Sampling ABC (RIS-ABC)**, directly applying the paper's concepts.
