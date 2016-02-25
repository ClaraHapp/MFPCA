#' The functional CP-TPA algorithm
#' 
#' This function implements the functional CP-TPA (FCP-TPA) algorithm, that 
#' calculates a smooth PCA for 3D tensor data (i.e. \code{N} observations of 2D 
#' images with dimension \code{S1 x S2}). The results are given in a 
#' CANDECOMP/PARAFRAC (CP) model format \deqn{X = sum_{k = 1}^K d_k u_k o v_k o 
#' w_k}  where \eqn{o} stands for the outer product. In this representation, the
#' outer product \eqn{v_k o w_k} can be regarded as the \eqn{k}-th eigenimage,
#' while \eqn{d_k u_k} represents the vector of individual scores for this
#' eigenimage and each observation.
#' 
#' The degree of smoothness is controlled by three smoothing parameters 
#' \eqn{\alpha_u, \alpha_v, \alpha_w}. There is hence one smoothing parameter 
#' for each direction, passed via the 
#' parameter \code{alphaVal}. 
#' 
#' @param X The data tensor of dimensions \code{N x S1 x S2}
#' @param K The number of eigentensors to be calculated
#' @param penMat A list of three roughness penalization matrices, one for each 
#'   direction (u, v & w).
#' @param alphaVal A list of three smoothing parameters along each 
#'   direction.
#' @param verbose Logical. If \code{TRUE} computational details are given on the
#'   standard output during calculation of the FPCA
#'  @param tol A numeric value, given the tolerance for relative and absolute values in the algorithm. Defaults to \code{1e-5}.
#'   
#' @return \item{d}{A vector of length \code{K}, containing the numeric weigths 
#'   \eqn{d_k} in the CP model} \item{U}{A matrix of dimensions \code{N x K}, 
#'   containing the eigenvectors \eqn{u_k} in the first dimension} \item{V}{A 
#'   matrix of dimensions \code{S1 x K}, containing the eigenvectors \eqn{v_k} 
#'   in the second dimension} \item{W}{A matrix of dimensions \code{S2 x K}, 
#'   containing the eigenvectors \eqn{w_k} in the third dimension.}
#'   
#' @references G. I. Allen, "Multi-way Functional Principal Components 
#'   Analysis", IEEE International Workshop on Computational Advances in 
#'   Multi-Sensor Adaptive Processing, 2013.
#'   
#' @seealso fcptpaBasis
FCP_TPA <- function(X, K, penMat, alphaVal, verbose = FALSE, tol = 1e-5)
{
  dimX <- dim(X)
  
  # initialize the norm-one vectors u,v,w randomly
  u <- runif(dimX[1], min = -1, max = 1); u <- u/normVec(u)
  v <- runif(dimX[2], min = -1, max = 1); v <- v/normVec(v)
  w <- runif(dimX[3], min = -1, max = 1); w <- w/normVec(w)
  
  # initialize smoothness parameters and smoothing matrices
  Iu <- diag(dimX[1])
  Iv <- diag(dimX[2])
  Iw <- diag(dimX[3])
  
  # initialize results
  d <- rep(NA, K)
  
  U <- matrix(NA, ncol = K, nrow = dimX[1])
  V <- matrix(NA, ncol = K, nrow = dimX[2])
  W <- matrix(NA, ncol = K, nrow = dimX[3])
  
  for(k in 1:K)
  {
    if(verbose)
      cat("\nk = ", k, "\n")
    
    # initialize "old" versions
    uOld <- 0*u
    vOld <- 0*v
    wOld <- 0*w
    
    # repeat until convergence
    while(any(c(normVec(u - uOld)/normVec(u), normVec(v - vOld)/normVec(v), normVec(w - wOld)/normVec(w)) > tol))
    {
      # update u
      uOld <- u
      u <- solve(Iu + alphaVal$u * penMat$u, ttv(X, list(v, w), dim = 2:3)) /
        as.numeric(crossprod(v, v + alphaVal$v*penMat$v %*%v )*crossprod(w, w + alphaVal$w * penMat$w %*% w))#as.numeric(crossprod(v, solve(Sv, v))* crossprod(w, solve(Sw, w)))
      
      # update v
      vOld <- v
      v <- solve(Iv + alphaVal$v * penMat$v, ttv(X, list(u, w), dim = c(1,3)))/
        as.numeric(crossprod(u, u + alphaVal$u * penMat$u %*% u)*crossprod(w, w + alphaVal$w * penMat$w %*% w))
         
      # update w
      wOld <- w
      w <- solve(Iw + alphaVal$w * penMat$w, ttv(X, list(u, v), dim = 1:2)) /
        as.numeric(crossprod(u,  u + alphaVal$u * penMat$u %*% u)* crossprod(v, v + alphaVal$v*penMat$v %*% v ))
    }
    
    if(verbose)
      cat("Absolute error:\n
          u: ", normVec(u - uOld), ", v: ", normVec(v - vOld), "w: ",  normVec(w - wOld), "\n")
    
    # scale vectors to have norm one
    u <- u/normVec(u)
    v <- v/normVec(v)
    w <- w/normVec(w)
    
    # calculcate results
    d[k] <- ttv(X, list(u,v,w), dim = 1:3)
    U[,k] <- u
    V[,k] <- v
    W[,k] <- w
    
    # update X
    X <- X - d[k] * u %o% v %o% w #u[,1] %o% v[,1] %o% w[,1] # outer product   
  }
  
  return(list(d = d, U = U, V = V, W = W))
}

#' Calculate the euclidean norm of a vector
#' 
#' @param x The vector for which the norm is to be calculated
#' 
#' @return The euclidean norm of \code{x}.
#' 
#' @keywords internal
normVec <- function(x) sqrt(sum(x^2))

#' Generalized cross-validation for the FCP-TPA algorithm
#' 
#' These functions calculate the generalized cross-validation criteria for the 
#' smoothing parameters \eqn{\alpha_v} and \eqn{\alpha_w} that are calculated within the 
#' \code{\link{FCP_TPA}} algorithm.
#' 
#' @section Warning: In test scenarios, the GCV always tended to select the 
#'   highest possible value for \eqn{\alpha_v, \alpha_w}, i.e. it seemed to be a
#'   monotonically decreasing function in \eqn{\alpha_v, \alpha_w}.
#'   
#' @param alphaV,alphaW The values of the smoothing parameters \eqn{\alpha_v}, 
#'   \eqn{\alpha_w} (numerical).
#' @param data The tensor containing the data, an array of dimensions \code{N x 
#'   S1 x S2}.
#' @param u,v,w The current value of the eigenvectors \eqn{u_k, v_k, w_k} (not 
#'   normalized) of dimensions \code{N}, \code{S1} and \code{S2}, respectively.
#' @param penMat A list of three roughness penalization matrices, one for each 
#'   direction (u, v & w). See \code{\link{FCP_TPA}}.
#'   
#' @return The value of the GCV criterion.
#' 
#' @references G. I. Allen, "Multi-way Functional Principal Components 
#'   Analysis", IEEE International Workshop on Computational Advances in 
#'   Multi-Sensor Adaptive Processing, 2013.
#'   
#' @keywords internal
#'   
#' @seealso FCP_TPA
gcvV <- function(alphaV, data, u, v, w, penMat, alphaU, alphaW)
{
  m <- length(v)
  
  res <- 1/m* normVec(ttv(data, list(u,w), c(1,3))/(normVec(u) * normVec(w))^2 - v)^2 / 
    (1 - (normVec(u) * normVec(w))^2/m * sum(1/eigen(diag(m) + alphaV*penMat$v)$values) / (as.numeric(crossprod(u, u + alphaU * penMat$u %*% u)) * as.numeric(crossprod(w, w + alphaW * penMat$w %*% w))) )^2  
  
  return(res)
}

#' @rdname gcvV
gcvW <- function(alphaW, data, u, v, w, penMat, alphaU, alphaV)
{
  p <- length(w)
  
  res <- 1/p* normVec(ttv(data, list(u,v), c(1,2))/(normVec(u) * normVec(v))^2 - w)^2 / 
    (1 - (normVec(u) * normVec(v))^2/p * sum(1 / eigen(diag(p) + alphaW * penMat$w)$values ) / ( as.numeric(crossprod(u, u + alphaU * penMat$u %*% u)) * as.numeric(crossprod(v, v + alphaV*penMat$v %*%v ))) )^2  
  
  return(res)
}
