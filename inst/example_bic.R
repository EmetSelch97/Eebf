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
