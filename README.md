# Exact Exploratory Bi-factor Analysis
This package provides the R code for the constraint-based exact exploratory bi-factor analysis method, including model estimation and model selection.

## How to install this package
To install this package, try remotes::install_github("EmetSelch97/Eebf").

## How to estimate a bi-factor structure when the number of factors is known
To estimate a bi-factor structure, refer to inst/example_estimation.R or 
```R
library(Eebf)
set.seed(2024)
J=15
G=3
n=2000

Q = Q_generator(J,G)

para_sim = parameter_generator(Q)
L_true = para_sim$L
Psi_true = para_sim$Psi
D_true = para_sim$D
Cov_true = para_sim$Cov

tlist = list(
  L_true = L_true,
  D_true = D_true,
  Psi_true = Psi_true
)


S=sampling_process(Cov_true,n)

results = Ebfa_solver(S=S,G=G)
Q_est = results$Q_est
L_est = results$L_est
Psi_est = results$Psi_est
D_est = results$D_est
nll_est = results$nll_est

EMC = Exact_matching(Q,Q_est)
ACC = Average_correct(Q,Q_est)
plist = list(
  L_est = L_est,
  D_est = D_est,
  Psi_est = Psi_est
)
MSE_results = MSE_parameter(plist,tlist)
```
## How to select the number of factors
To select the number of factors, refer to inst/example_bic.R or
```R
library(Eebf)
set.seed(2024)
J=15
G=3
n=2000

Q = Q_generator(J,G)

para_sim = parameter_generator(Q)
L_true = para_sim$L
Psi_true = para_sim$Psi
D_true = para_sim$D
Cov_true = para_sim$Cov

S=sampling_process(Cov_true,n)

results = bic_process((G-1),(G+1),S,n)
G_est = results$G_est
```
## Further details of this package
This package involves six R scripts.  
In utils.R, some basic functions used in the optimization step of this mehtod are provided.  
In opt_step.R, the model estimation method via the Augmented Lagrangian method is provided in function Ebfa_solver().  
In bic_step.R, the model section process via BIC is provided. By function bic_process(), users can select the number of factors of the bi-factor model.  
In generator.R, the data generating process for simulation studies are provided. Use Q_generator() to generate the bi-factor structure, parameter_generator() to generate bi-factor parameters including loading matrix, correlation matrix and unique variance matrix 
and sampling_process to generate sample covariance matrix.  
In eval.R, the evaluation criteria for the simulation studies are provided.  
