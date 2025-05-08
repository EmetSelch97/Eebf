#' The Augmented Lagrangian method for exploring an exact bi-factor model
#' @param x_init initial value for parameters
#' @param J the number of items
#' @param G the number of group factors
#' @param S the sample covariance matrix
#' @param Pair the constraint for bi-factor model
#' @param gamma_0 the initial value for gamma, a matrix-like ALM parameter
#' @param rho_0 the initial value for rho, a scalar, ALM parameter
#' @param rho_sigma a hyper-parameter for ALM method
#' @param theta a hyper-parameter for ALM method
#' @param tol the tolerance for ALM method
#' @param max_iter the maximal iteration number of the ALM method
#' @return A list containing:
#'   \itemize{
#'     \item \code{x_est}: the estimated parameters of the bi-factor model
#'     \item \code{iter_num}: the iteration number of the bi-factor model
#'     \item \code{L_est}: the estimated loading matrix
#'     \item \code{Psi_est}: the estimated correlation matrix
#'     \item \code{D_est}: the estimated unique variance matrix
#'     \item \code{Q_est}: the estimated bi-factor structure
#'   }
#' @export
#' @importFrom stats optim
alm_solve<-function(J,G,S,Pair,gamma_0,rho_0,rho_sigma,theta,tol,max_iter,x_init = numeric(0)){

  if (length(x_init) == 0 ){
    x_init = init_alm(J,G)
  }

  ep = J*(G+2) + G*(G-1)/2

  x_old = x_init
  gamma = gamma_0
  rho = rho_0
  old_dec = para_decompose(x_old,J,G)
  L_old = old_dec$L
  cons_val_old = cons_fun(L_old,J,G,Pair)

  dist_para = sqrt(sum(x_old^2))/sqrt(ep)
  dist_val = bf_dist(L_old)
  iter_num = 0
  while (max(dist_para,dist_val)>tol & iter_num<max_iter){
    result = optim(
      par = x_old,
      fn = objective_function,
      gr = alm_gd,
      method = "L-BFGS-B",
      J = J,
      G = G,
      S = S,
      gamma = gamma,
      rho = rho,
      Pair = Pair

    )
    x_new = result$par

    dist_para = sqrt(sum((x_old-x_new)^2))/sqrt(ep)

    new_dec = para_decompose(x_old,J,G)
    L_new = new_dec$L
    cons_val_new = cons_fun(L_new,J,G,Pair)

    gamma = gamma + rho*cons_val_new
    if (sqrt(sum(cons_val_new^2))>theta*sqrt(sum(cons_val_old^2))){
      rho = rho * rho_sigma
    }

    x_old = x_new
    old_dec = para_decompose(x_old,J,G)
    L_old = old_dec$L
    cons_val_old = cons_fun(L_old,J,G,Pair)
    dist_val = bf_dist(L_old)
    iter_num = iter_num +1
  }

  D_old = old_dec$D
  Psi_old = old_dec$Psi
  Q_est = bf_compress(L_old,tol)
  return(list(
    x_est = x_old,
    iter_num = iter_num,
    L_est = L_old,
    Psi_est = Psi_old,
    D_est = D_old,
    Q_est = Q_est
  ))
}

#' EBFA solver
#' @param S the sample covariance matrix
#' @param G the number of group factors
#' @param rho_0 the initial value of rho in ALM method, default is 1
#' @param rho_sigma the hyper-parameter in ALM method, default is 10
#' @param theta the hyper-parameter in ALM method, default is 0.25
#' @param tol the tolerance for the ALM method, default is 0.01
#' @param max_iter the maximal iteration number, default is 100
#' @param random_starts the number of random starting points, default is 50
#' @return A list containing:
#'   \itemize{
#'     \item \code{nll_est}: the negative-log-likelihood of the estimated parameters of the bi-factor model
#'     \item \code{L_est}: the estimated loading matrix
#'     \item \code{D_est}: the estimated unique variance matrix
#'     \item \code{Psi_est}: the estimated correlation matrix
#'     \item \code{Q_est}: the estimated bi-factor structure
#'   }
#' @export
Ebfa_solver<-function(S,G,rho_0=1,rho_sigma=10,theta=0.25,tol=0.01,max_iter=100,random_starts=50){
  J = nrow(S)
  Pair = Bf_cons_pair(G)
  p = length(Pair)
  gamma_0 = matrix(0,nrow = p,ncol=J)

  L_list = list()
  D_list = list()
  Psi_list = list()
  nll_list = list()
  for (i in 1:random_starts){
    x_init = init_alm(J,G)
    results = alm_solve(J,G,S,Pair,gamma_0,rho_0,rho_sigma,theta,tol,max_iter,x_init)
    if (results$iter_num < max_iter){
      L_alm = results$L_est
      Psi_alm = results$Psi_est
      D_alm = results$D_est
      nll_alm = nll(L_alm,Psi_alm,D_alm,S)

      L_list = c(L_list,list(L_alm))
      D_list = c(D_list,list(D_alm))
      Psi_list = c(Psi_list,list(Psi_alm))
      nll_list = c(nll_list,list(nll_alm))
    }
  }
  restart = 0
  while (length(nll_list)<(random_starts/2)){
    restart = restart + 1
    for (i in 1:random_starts){
      x_init = init_alm(J,G)
      results = alm_solve(J,G,S,Pair,gamma_0,rho_0,rho_sigma,theta,tol,max_iter,x_init)
      if (results$iter_num < max_iter){
        L_alm = results$L_est
        Psi_alm = results$Psi_est
        D_alm = results$D_est
        nll_alm = nll(L_alm,Psi_alm,D_alm,S)

        L_list = c(L_list,list(L_alm))
        D_list = c(D_list,list(D_alm))
        Psi_list = c(Psi_list,list(Psi_alm))
        nll_list = c(nll_list,list(nll_alm))
      }
    }
  }
  if (length(nll_list)>0){
    nll_vec = unlist(nll_list)
    min_post = which.min(nll_vec)

    nll_est =  nll_list[[min_post]]
    L_est = L_list[[min_post]]
    D_est = D_list[[min_post]]
    Psi_est = Psi_list[[min_post]]
    Q_est = bf_compress(L_est,tol)
  } else{
    stop('No converged solutions')
  }


  return(list(
    nll_est = nll_est,
    L_est = L_est,
    D_est = D_est,
    Psi_est = Psi_est,
    Q_est = Q_est
  ))
}
