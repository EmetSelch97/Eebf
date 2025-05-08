#' negative-log-likelihood function
#'
#' @param L input of loading matrix
#' @param Phi input of correlation matrix
#' @param D input of diagonal matrix
#' @param S sample covariance matrix
#' @return negative-log-likelihood of data
#' @export

nll<-function(L,Phi,D,S){
  Cov = L %*% (Phi %*% t(L)) +D
  R = solve(Cov,S)
  loss = log(2*pi*det(Cov))/2 + sum(diag(R))/2
  return(loss)
}

#' re-parameterization of Phi
#'
#' @param x input of free parameters in the correlation matrix
#' @param G the number of group factors
#' @return the Cholesky decomposition of the correlation matrix
#' @export

stan_trans<-function(x,G){
  if (length(x) != G*(G-1)/2){
    stop('The dimension of input does not match the number of free parameters in Phi')
  }
  h = matrix(0,nrow = G,ncol = G)
  for (i in 2:G){
    h[1:(i-1),i] = x[(1+(i-1)*(i-2)/2):((i-1)*i/2)]
  }
  Z = tanh(h)
  U = matrix(0,nrow=G,ncol=G)
  U[1,1] = 1
  U[1,2:G] = Z[1,2:G]
  for (i in 2:G){
    U[i,i] = U[(i-1),i]*sqrt(1-Z[(i-1),i]^2)/Z[(i-1),i]
    if (i<G){
      U[i,(i+1):G] = Z[i,(i+1):G]*U[(i-1),(i+1):G]*sqrt(1-Z[(i-1),(i+1):G]^2)/Z[(i-1),(i+1):G]
    }

  }
  return(U)
}

#' transformation of free parameters into the form of loading matrix, correlation matrix and unique variance matrix
#' @param x input of free parameters of the full bi-factor model
#' @param J the number of items
#' @param G the number of group factors
#' @return A list containing:
#'   \itemize{
#'     \item \code{L}: the loading matrix
#'     \item \code{D}: the unique variance matrix
#'     \item \code{d}: the parameter of input satisfying diag(d^2) = D
#'     \item \code{Phi}: the Cholesky decomposition of the correlation matrix
#'     \item \code{Psi}: the correlation matrix
#'     \item \code{Cov}: the true covariance matrix
#'   }
#' @export
para_decompose<-function(x,J,G){
  if (length(x) != (J*(2+G) + G*(G-1)/2)){
    stop('The dimension of input x does not match the number of free parameters of the bi-factor model')
  }
  Phi = matrix(0,nrow = 1+G, ncol = 1+G)
  Phi[1,1] = 1
  Phi[2:(1+G),2:(1+G)] = stan_trans(x[1:(G*(G-1)/2)],G)

  Psi = t(Phi) %*% Phi

  L = matrix(x[(1+G*(G-1)/2): (J*(1+G)+G*(G-1)/2)],nrow = J,ncol=(G+1))
  d = x[(1+J*(1+G)+G*(G-1)/2):(J*(2+G) + G*(G-1)/2)]
  D = diag(d^2)

  Cov = L %*% (Psi %*% t(L)) + D

  return(list(
    L = L,
    D = D,
    d = d,
    Phi = Phi,
    Psi = Psi,
    Cov = Cov
  ))
}

#' chain rule of the gradient of Phi
#' @param x the input vector
#' @param dU the gradient in the chain rule
#' @param G the number of group factors
#' @return a intermediate result for calculating the gradient of the ALM method
#' @export
gd_rp <-function(x,dU,G){
  if (length(x) != G * (G - 1) / 2) {
    stop("The dimension of x does not match the number of free parameters of Phi")
  }
  if (length(dU) != (G - 1) * (G + 2) / 2) {
    stop("The dimension of dU does not match the dimension of matrix in the chain rule")
  }
  Z <- tanh(x)
  dZ <- 1 - Z^2

  n_row <- G * (G - 1) / 2
  n_col <- (G - 1) * (G + 2) / 2
  A <- matrix(0, nrow = n_row, ncol = n_col)

  for (i in 1:(G-1)){
    y = Z[(1+i*(i-1)/2): (i*(i+1)/2)]
    ly = length(y)
    sy = sqrt(1-y^2)
    inv_sy = 1 / (sy + 1e-15)

    scaler_full = rep(1, ly)
    if (ly>1) {
      scaler_full[-1] = sy[-ly]
    }
    scaler_full <- cumprod(scaler_full)

    A_sub <- matrix(0, nrow = i, ncol = i + 1)

    for (j in 1:i) {
      A_sub[j, j] = 1

      rl = i - j + 1
      r_1 = rep(1, rl)
      r_2 = rep(1, rl)
      if (rl > 1) {
        r_1[-rl] = y[(j+1):ly]

        r_2[-1] = sy[(j+1):ly]
        r_2 = cumprod(r_2)
      }

      product = -r_1 * r_2 * inv_sy[j] * y[j]

      if (length(product) > 0) {
        A_sub[j, (j + 1):(i + 1)] = product
      }
    }

    scaler_diag = diag(scaler_full)
    A_sub_scaled = scaler_diag %*% A_sub

    A[(1+i*(i-1)/2):(i*(i+1)/2), (i*(i+1)/2):(i*(i+3)/2)] = A_sub_scaled
  }
  gd <- (A %*% dU) * dZ

  return(gd)
}

#' initial value of bi-factor optimization
#' @param J the number of items
#' @param G the number of group factors
#' @return the initial value of bi-factor optimization
#' @export
#' @importFrom stats runif
init_alm<-function(J,G){
  tl = J*(2+G) + G*(G-1)/2
  x = numeric(tl)
  l1 = G*(G-1)/2
  l2 = G*(G-1)/2 + J
  l3 = J*(1+G) + G*(G-1)/2
  x[1:l1] = runif(l1, min = -0.5, max = 0.5)
  x[(1+l1):l2] = runif(J)
  x[(1+l2):l3] = runif(G * J, min = -0.5, max = 0.5)
  x[(1+l3):tl] = 1+runif(J)
  return(x)
}

#' bi-factor approximation evaluation
#' @param L the input loading matrix
#' @return the distance of the given loading matrix to the space of bi-factor loading matrix
#' @export
bf_dist<-function(L){
  G = ncol(L)-1
  sorted_L = t(apply(abs(L[,2:(G+1)]), 1, sort, decreasing = TRUE))
  dist = max(sorted_L[,2])
  return(dist)
}

#' compress a bi-factor structure from loading
#' @param L the input loading matrix
#' @param tol the tolerance of the absolute value of the loading matrix
#' @return a 0-1 matrix, indicating the bi-factor structure learned from the loading matrix
#' @export
bf_compress<-function(L,tol){
  r = ncol(L)
  J = nrow(L)
  Q = matrix(0,nrow=J,ncol=r)
  Q[,1] = rep(1,J)
  Q[,2:r] = as.integer(abs(L[,2:r])>tol)
  return(Q)
}
