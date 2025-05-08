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




