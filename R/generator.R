#' generate a bi-factor structure as an example
#' @param J the number of item
#' @param G the number of factors
#' @return a 0-1 matrix indicating the bi-factor structure
#' @export
Q_generator<-function(J,G){
  if (J%%G == 0){
    p = J%/%G
    Q = matrix(0,nrow=J,ncol=(G+1))
    Q[,1] = replicate(J,1)
    for (i in 1:p){
      Q[(1+(i-1)*G):(i*G),2:(G+1)] = diag(G)
    }
  } else{
    stop('The case of (J,G) is not considered in the example')
  }
  return(Q)
}

#' generate the loading matrix, correlation matrix ,unique variance matrix and covariance matrix given factor structure
#' @param Q the factor structure matrix
#' @param case the case for generating the parameters, 'exact' indicates an exact bi-factor structure and 'approx' indicates an approximate bi-factor structure
#' @return A list containing:
#'   \itemize{
#'     \item \code{L}: the true loading matrix
#'     \item \code{Psi}: the true correlation matrix
#'     \item \code{D}: the true unique variance matrix
#'     \item \code{Cov}: the true covariance matrix
#'   }
#' @export
#' @importFrom stats runif rbinom
parameter_generator<-function(Q,case='exact'){
  J = nrow(Q)
  C = ncol(Q)
  G = C-1
  if (case == 'exact'){
    D = diag(J)
    Bin = matrix(1 - 2 * rbinom(J * G, size = 1, prob = 0.5), nrow = J, ncol = G)
    Sign = matrix(1, nrow = J, ncol = C)
    Sign[, 2:(G+1)] = Bin

    L1 = matrix(runif(J * (1 + G), min = 0.2, max = 1), nrow = J, ncol = C)
    L = L1*Sign*Q

    Phi = matrix(0,nrow=C,ncol = C)
    Phi[1,1] = 1
    phi = runif(G*(G-1)/2,min=-0.5,max=0.5)
    Phi[2:C,2:C] = stan_trans(phi,G)
    Psi = t(Phi) %*% Phi

    Cov = L %*% (Psi %*% t(L)) + D
  } else if (case == 'approx'){
    D = diag(J)
    Bin = matrix(1 - 2 * rbinom(J * G, size = 1, prob = 0.5), nrow = J, ncol = G)
    Sign = matrix(1, nrow = J, ncol = C)
    Sign[, 2:(G+1)] <- Bin

    L1 = matrix(runif(J * (1 + G), min = 0.2, max = 1), nrow = J, ncol = C)
    L2 = matrix(runif(J * (1 + G), min = 0, max = 0.1), nrow = J, ncol = C)
    L = ifelse(Q==1,L1*Sign,L2*Sign)

    Phi = matrix(0,nrow=C,ncol = C)
    Phi[1,1] = 1
    phi = runif(G*(G-1)/2,min=-0.5,max=0.5)
    Phi[2:C,2:C] = stan_trans(phi,G)
    Psi = t(Phi) %*% Phi

    Cov = L %*% (Psi %*% t(L)) + D
  } else{
    stop('The case is not considered in the examples')
  }
  return(list(
    L = L,
    Psi = Psi,
    D = D,
    Cov = Cov
  ))
}


#' generate the sample covariance matrix given the covariance matrix
#' @param Cov the covariance matrix
#' @param n the sample size
#' @return the sample covariance matrix
#' @export
#' @importFrom MASS mvrnorm
sampling_process<-function(Cov,n){
  J = nrow(Cov)
  mu = numeric(J)
  samples = mvrnorm(n=n,mu=mu,Sigma = Cov)
  S = t(samples) %*% samples /n
  return(S)
}
