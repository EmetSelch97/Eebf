# Selecting the number of factors by Bayesian information criterion
#'
#' @title Model Selection via BIC
#' @param G_min the lower bound for the number of factors
#' @param G_max the upper bound for the number of factors
#' @param S the sample covariance matrix
#' @param n the sample size
#' @param rho_0 the initial value of rho in ALM method, default is 1
#' @param rho_sigma the hyper-parameter in ALM method, default is 10
#' @param theta the hyper-parameter in ALM method, default is 0.25
#' @param tol the tolerance for the ALM method, default is 0.01
#' @param max_iter the maximal iteration number, default is 100
#' @param random_starts the number of random starting points, default is 50
#' @return A list containing:
#'   \itemize{
#'     \item \code{G_est}: the selected number of group factors based on BIC
#'     \item \code{L_est}: the estimated loading matrix given G_est group factors
#'     \item \code{D_est}: the estimated unique variance matrix given G_est group factors
#'     \item \code{Psi_est}: the estimated correlation matrix given G_est group factors
#'     \item \code{Q_est}: the estimated bi-factor structure given G_est group factors
#'   }
#' @export
bic_process<-function(G_min,G_max,S,n,rho_0=1,rho_sigma=10,theta=0.25,tol=0.01,max_iter=100,random_starts=50){
  Len = G_max-G_min+1
  BIC_list = numeric(Len)
  L_list = vector("list",Len)
  Psi_list = vector("list",Len)
  D_list = vector("list",Len)
  for (g in G_min : G_max){
    results = Ebfa_solver(S,g,rho_0,rho_sigma,theta,tol,max_iter,random_starts)
    nll_bic = results$nll_est
    L_bic = results$L_est
    D_bic = results$D_est
    Psi_bic = results$Psi_est

    Q_bic = bf_compress(L_bic,tol)

    bic = 2*n*nll_bic + log(n)*(g-1)*g/2
    BIC_list[g-G_min+1] = bic

    L_list[[g-G_min+1]] = L_bic
    Psi_list[[g-G_min+1]] = Psi_bic
    D_list[[g-G_min+1]] = D_bic
  }
  best_post = which.min(BIC_list)
  g_est = best_post + G_min -1
  L_est = L_list[[best_post]]
  D_est = D_list[[best_post]]
  Psi_est = Psi_list[[best_post]]

  Q_est = bf_compress(L_est,tol)
  return(list(
    G_est = g_est,
    L_est = L_est,
    D_est = D_est,
    Psi_est = Psi_est,
    Q_est = Q_est
  ))
}
