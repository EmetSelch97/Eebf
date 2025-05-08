#' Bi-factor permutation
#' @param G the number of group factors of a bi-factor model
#' @return list of all possible permutation matrix of the group factors
#' @export
#' @importFrom gtools permutations
BF_permutation <- function(G) {
  elements <- 1:G
  all_perms <- permutations(G, G, v = elements)

  n_row = nrow(all_perms)
  Mat_list = vector('list',n_row)
  for (i in 1:n_row){
    mat = matrix(0,nrow=G,ncol = G)
    for (j in 1:G){
      mat[j,all_perms[i,j]] = 1
    }
    Mat_list[[i]] = mat
  }

  return(Mat_list)
}

#' Exact Match Criterion for bi-factor structure
#'
#' @param Q the default bi-factor structure
#' @param Q_est the estimated bi-factor structure
#' @param Rot_list the rotation matrix considered in the criterion, default is an empty list
#' @return a 0-1 variable, which equals to 1 when the bi-factor structure is correctly learned and 0 otherwise
#' @export
Exact_matching<-function(Q,Q_est,Rot_list = list()){
  if (nrow(Q)!=nrow(Q_est) | ncol(Q)!=ncol(Q_est)){
    stop('The dimensions of Q_est and Q do not match')
  } else{
    G = ncol(Q)-1
  }
  if (length(Rot_list) ==0){
    Rot_list = BF_permutation(G)
  }
  p = length(Rot_list)
  EMC_list = 1+numeric(p)
  for (i in 1:p){
    Q_rot = Q
    Q_rot[,2:(G+1)] = Q_rot[,2:(G+1)] %*% Rot_list[[i]]
    EMC_list[i] = sum(Q_rot != Q_est)
  }
  EMC = 0
  if (min(EMC_list) ==0){
    EMC=1
  }
  return(EMC)
}

#' Calculate the average intersection of two bi-factor structure
#'
#' @param Q the default bi-factor structure
#' @param Q_est the estimated bi-factor structure
#' @return the average intersection of two bi-factor structure ranges between 0 to 1
#' @export
Compare_Q<-function(Q,Q_est){
  J = nrow(Q)
  r = ncol(Q)
  if (nrow(Q_est) != J | ncol(Q_est) != r){
    stop('The dimension of Q_est does not match the dimension of Q')
  }

  average_cr = 0
  for (i in 2:r){
    zero_post = which(Q[,i] ==0)
    zero_est_post = which(Q_est[,i] ==0)
    zero_inter = intersect(zero_post,zero_est_post)

    one_post = which(Q[,i] ==1)
    one_est_post = which(Q_est[,i] ==1)
    one_inter = intersect(one_post,one_est_post)

    average_cr= average_cr + length(zero_inter) +length(one_inter)
  }
  return(average_cr/(J*(r-1)))
}

#' Average Correctness Criterion for bi-factor structure
#'
#' @param Q the default bi-factor structure
#' @param Q_est the estimated bi-factor structure
#' @param Rot_list the rotation matrix considered in the criterion, default is an empty list
#' @return value ranging between 0-1, a larger value indicates a higher correctness
#' @export
Average_correct<-function(Q,Q_est,Rot_list=list()){
  G = ncol(Q)-1
  if (length(Rot_list) ==0){
    Rot_list = BF_permutation(G)
  }
  p = length(Rot_list)
  ACC_list = 1+numeric(p)

  for (i in 1:p){
    Q_rot = Q
    Q_rot[,2:(1+G)] = Q_rot[,2:(1+G)] %*% Rot_list[[i]]
    ACC_list[i] = Compare_Q(Q_rot,Q_est)
  }
  return(max(ACC_list))
}



#' Calculate the mean squared error of bi-factor parameters after rotations
#' @param plist A list containing:
#'   \itemize{
#'     \item \code{L_est}: the estimated loading matrix
#'     \item \code{D_est}: the estimated unique variance matrix
#'     \item \code{Psi_est}: the estimated correlation matrix
#'   }
#' @param tlist A list containing:
#'   \itemize{
#'     \item \code{L_true}: the true loading matrix
#'     \item \code{D_true}: the true unique variance matrix
#'     \item \code{Psi_true}: the true correlation matrix
#'   }
#' @param Rot_list the rotation list for evaluation
#' @return A list containing:
#'   \itemize{
#'     \item \code{MSE_L}: the MSE of loading matrix
#'     \item \code{MSE_D}: the MSE of unique variance matrix
#'     \item \code{MSE_Psi}: the MSE of correlation matrix
#'   }
#' @export
MSE_parameter<-function(plist,tlist,Rot_list = list()){
  L_est = plist$L_est
  D_est = plist$D_est
  Psi_est = plist$Psi_est

  L_true = tlist$L_true
  D_true = tlist$D_true
  Psi_true = tlist$Psi_true

  G_est = ncol(L_est)-1
  G = ncol(L_true)-1
  J = nrow(L_est)

  if (G_est != G){
    stop('The dimensions of loading matrix do not match')
  }

  if (length(Rot_list) == 0){
    Rot_list = BF_permutation(G)
  }

  MSE_D = sum((D_est-D_true)^{2})/J
  p = length(Rot_list)

  L_list = 1+numeric(p)

  Psi_list = 1+numeric(p)
  for (i in 1:p){
    Rot = matrix(0,nrow = (G+1),ncol = (G+1))
    Rot[1,1] = 1
    Rot[2:(G+1),2:(G+1)] = Rot_list[[i]]

    L_rot = L_true %*% Rot
    Sign = diag(diag(sign(t(L_est) %*% L_rot)))
    L_list[i] = sum((L_est - L_rot %*% Sign)^2)/(J*(1+G))

    Psi_rot = Sign %*% (t(Rot) %*% (Psi_true %*% (Rot %*% Sign)))
    Psi_list[i] = sum((Psi_est - Psi_rot)^2)/((1+G)^2)
  }
  MSE_L = min(L_list)
  Loc = which.min(L_list)
  MSE_Psi = Psi_list[Loc]
  return(list(
    MSE_L = MSE_L,
    MSE_D = MSE_D,
    MSE_Psi = MSE_Psi
  ))
}
