#' bi-factor constraint pair of columns
#'
#' @param G the number of group factors
#' @return A list containing all pairs of columns in the equality constraints
#' @export
#' @importFrom utils combn

Bf_cons_pair <- function(G) {
  combn(2:(G+1), 2, simplify = FALSE)
}

#' bi-factor constraint value function
#'
#' @param L the input loading matrix
#' @param J the number of items
#' @param G the number of group factors
#' @param Pair the list indicating which pair of columns have constraints
#' @return a matrix consisting of the values of equality constraints
#' @export

cons_fun <-function(L,J,G,Pair){
  p = length(Pair)
  ep = ((G-1)*G)/2
  if (p !=ep){
    stop("The dimension of constraint pair does not match bi-factor model")
  }
  cons_val <- matrix(0, nrow = p, ncol = J)
  for (i in 1:p){
    Z = Pair[[i]]
    a = Z[1]
    b = Z[2]
    cons_val[i,] = L[,a]*L[,b]
  }
  return(cons_val)
}

#' bi-factor constraint gradient function
#'
#' @param L the loading matrix
#' @param J the number of items
#' @param G the number of group factors
#' @param gamma a matrix-like ALM parameter
#' @param rho a scalar-like ALM parameter
#' @param Pair a list containing all pairs of columns in the equality constraints
#' @return a vector, which is the gradient of the bi-factor constraint value function
#' @export
cons_gd <-function(L,J,G,gamma,rho,Pair){
  p = length(Pair)
  ep = (G-1)*G/2
  if (p !=ep){
    stop("The dimension of constraint pair does not match bi-factor model")
  }
  L2 = L^2
  grad_L <- matrix(0, nrow = J, ncol = G + 1)
  for (i in 1:ep){
    Z = Pair[[i]]
    a = Z[1]
    b = Z[2]
    gm = as.vector(gamma[i, ])

    grad_L[, a] = grad_L[, a] + gm * L[, b] + rho * L[, a] * L2[, b]
    grad_L[, b] <- grad_L[, b] + gm * L[, a] + rho * L[, b] * L2[, a]
  }
  grad = numeric(ep+(G+2)*J)
  grad[(1+ep):(ep+(G+1)*J)] = as.vector(grad_L)
  return(grad)
}

#' Objective function in the Augmented Lagrangian method
#' @param x the input free parameters in bi-factor model
#' @param J the number of items
#' @param G the number of group factors
#' @param S the sample covariance matrix
#' @param gamma a matrix-like ALM parameter
#' @param rho a scalar-like ALM parameter
#' @param Pair a list containing all pairs of columns in the equality constraints
#' @return the objective function value in the Augmented Lagrangian method
#' @export
objective_function<-function(x,J,G,S,gamma,rho,Pair){

  ep = G*(G-1)/2 + (G+2)*J
  if (length(x) != ep){
    stop('The dimension of input x does not match the number of free parameters in bi-factor model in the objective function')
  }
  x_dc = para_decompose(x,J,G)

  L = x_dc$L
  D = x_dc$D
  Psi = x_dc$Psi

  loss_part = nll(L,Psi,D,S)

  cons_val = cons_fun(L,J,G,Pair)
  pen_part = sum(gamma*cons_val) + 0.5*rho*sum(cons_val^2)

  return(loss_part+pen_part)
}

#' gradient function in the Augmented Lagrangian method
#' @param x input free parameters in bi-factor model
#' @param J the number of items
#' @param G the number of group factors
#' @param S the sample covariance matrix
#' @param gamma a matrix-like ALM parameter
#' @param rho a scalar-like ALM parameter
#' @param Pair a list containing all pairs of columns in the equality constraints
#' @return the gradient of the objective function in the Augmented Lagrangian method
#' @export
alm_gd<-function(x,J,G,S,gamma,rho,Pair){

  ep = G*(G-1)/2 + (G+2)*J
  if (length(x) != ep){
    stop('The dimension of input x does not match the number of free parameters in bi-factor model in the gradient function')
  }

  x_dc = para_decompose(x,J,G)

  L = x_dc$L
  D = x_dc$D
  d = x_dc$d
  Phi = x_dc$Phi
  Psi = x_dc$Psi
  Cov = x_dc$Cov

  cons_grad = cons_gd(L,J,G,gamma,rho,Pair)

  nll_grad = numeric(ep)

  inv_Cov = solve(Cov)

  Sdw = inv_Cov %*% (S %*% inv_Cov)

  diff = inv_Cov-Sdw

  LP = L %*% Psi
  PL = Phi %*% t(L)

  dL = diff %*% LP
  nll_grad[(1+G*(G-1)/2):(G*(G-1)/2+(G+1)*J)] = as.vector(dL)

  dd = diag(diff) * d
  nll_grad[(1+G*(G-1)/2+(G+1)*J):ep] = dd

  dPhi = PL %*% (diff %*% L)
  dU = numeric((G-1)*(G+2)/2)
  for (i in 1:(G-1)){
    dU[(i*(i+1)/2):(i*(i+3)/2)] = dPhi[2:(i+2),i+2]
  }
  dphi = gd_rp(x[1:(G*(G-1)/2)],dU,G)
  nll_grad[1:(G*(G-1)/2)] = dphi
  return(nll_grad + cons_grad)
}
