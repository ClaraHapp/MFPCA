#' The functional CP-TPA algorithm
#' 
#' This function implements the functional CP-TPA (FCP-TPA) algorithm, that 
#' calculates a smooth PCA for 3D tensor data (i.e. \code{N} observations of 2D 
#' images with dimension \code{S1 x S2}). The results are given in a 
#' CANDECOMP/PARAFRAC (CP) model format \deqn{X = \sum_{k = 1}^K d_k \cdot u_k 
#' \circ v_k \circ w_k}{X = \sum d_k u_k \%o\% v_k \%o\% w_k}  where 
#' \eqn{\circ}{\%o\%} stands for the outer product, \eqn{d_k} is a scalar and 
#' \eqn{u_k, v_k, w_k} are eigenvectors for each direction of the tensor. In 
#' this representation, the outer product \eqn{v_k \circ w_k}{v_k \%o\% w_k} can
#' be regarded as the \eqn{k}-th eigenimage, while \eqn{d_k \cdot u_k}{d_k u_k} 
#' represents the vector of individual scores for this eigenimage and each 
#' observation.
#' 
#' The smoothness of the eigenvectors \eqn{v_k, w_k} is induced by penalty 
#' matrices for both image directions, that are weighted by smoothing parameters
#' \eqn{\alpha_{vk}, \alpha_{wk}}. The eigenvectors \eqn{u_k} are not smoothed, 
#' hence the algorithm does not induce smoothness along observations.
#' 
#' Optimal smoothing parameters are found via a nested generalized cross 
#' validation. In each iteration of the TPA (tensor power algorithm), the GCV 
#' criterion is optimized via \code{\link[stats]{optimize}} on the interval 
#' specified via \code{alphaRange$v} (or \code{alphaRange$w}, respectively).
#' 
#' The FCP_TPA algorithm is an iterative algorithm. Convergence is assumed if 
#' the relative difference between the actual and the previous values are all 
#' below the tolerance level \code{tol}. The tolerance level is increased 
#' automatically, if the algorithm has not converged after \code{maxIter} steps 
#' and if \code{adaptTol = TRUE}. If the algorithm did not converge after
#' \code{maxIter} steps (or \code{2 * maxIter}) steps, the function throws a
#' warning.
#' 
#' @param X The data tensor of dimensions \code{N x S1 x S2}.
#' @param K The number of eigentensors to be calculated.
#' @param penMat A list with entries \code{v} and \code{w}, containing a 
#'   roughness penalty matrix for each direction of the image. The algorithm 
#'   does not induce smoothness along observations (see Details).
#' @param alphaRange A list of length 2 with entries \code{v} and \code{w} , 
#'   containing the range of smoothness parameters to test for each direction.
#' @param verbose Logical. If \code{TRUE}, computational details are given on 
#'   the standard output during calculation of the FCP_TPA.
#' @param tol A numeric value, giving the tolerance for relative error values in
#'   the algorithm. Defaults to \code{1e-4}. It is automatically multiplied by 
#'   10 after \code{maxIter} steps, if \code{adaptTol = TRUE}.
#' @param maxIter A numeric value, the maximal iteration steps. Can be doubled, 
#'   if \code{adaptTol = TRUE}.
#' @param adaptTol Logical. If \code{TRUE}, the tolerance is adapted (multiplied
#'   by 10), if the algorithm has not converged after \code{maxIter} steps and 
#'   another \code{maxIter} steps are allowed with the increased tolerance, see 
#'   Details. Use with caution. Defaults to \code{TRUE}.
#'   
#' @return \item{d}{A vector of length \code{K}, containing the numeric weights 
#'   \eqn{d_k} in the CP model.} \item{U}{A matrix of dimensions \code{N x K}, 
#'   containing the eigenvectors \eqn{u_k} in the first dimension.} \item{V}{A 
#'   matrix of dimensions \code{S1 x K}, containing the eigenvectors \eqn{v_k} 
#'   in the second dimension.} \item{W}{A matrix of dimensions \code{S2 x K}, 
#'   containing the eigenvectors \eqn{w_k} in the third dimension.}
#'   
#' @references G. I. Allen, "Multi-way Functional Principal Components 
#'   Analysis", IEEE International Workshop on Computational Advances in 
#'   Multi-Sensor Adaptive Processing, 2013.
#'   
#' @seealso \code{\link{fcptpaBasis}}
#' 
#' @importFrom stats runif
#' 
#' @export FCP_TPA
#' 
#' @examples
#'  # set.seed(1234)
#'  
#'  N <- 100
#'  S1 <- 75
#'  S2 <- 75
#' 
#'  # define "true" components
#'  v <- sin(seq(-pi, pi, length.out = S1))
#'  w <- exp(seq(-0.5, 1, length.out = S2))
#'  
#'  # simulate tensor data with dimensions N x S1 x S2
#'  X <- rnorm(N, sd = 0.5) %o% v %o% w
#'  
#'  # create penalty matrices (penalize first differences for each dimension)
#'  Pv <- crossprod(diff(diag(S1)))
#'  Pw <- crossprod(diff(diag(S2)))
#'  
#'  # estimate one eigentensor
#'  res <- FCP_TPA(X, K = 1, penMat = list(v = Pv, w = Pw),
#'              alphaRange = list(v = c(1e-4, 1e4), w = c(1e-4, 1e4)),
#'              verbose = TRUE)
#'  
#'  # plot the results and compare to true values
#'  plot(res$V)
#'  points(v/sqrt(sum(v^2)), pch = 20)
#'  legend("topleft", legend = c("True", "Estimated"), pch = c(20, 1))
#'  
#'  plot(res$W)
#'  points(w/sqrt(sum(w^2)), pch = 20)
#'  legend("topleft", legend = c("True", "Estimated"), pch = c(20, 1))
FCP_TPA <- function(X, K, penMat, alphaRange, verbose = FALSE, tol = 1e-4, maxIter = 15, adaptTol = TRUE)
{
  ### check input parameters
  check_FCP_TPA_input(X, K, penMat, alphaRange, verbose, tol, maxIter, adaptTol)
  
  dimX <- dim(X)
  
  ### initialize all relevant values
  
  # initialize the norm-one vectors u, v, w randomly
  u <- stats::runif(dimX[1], min = -1, max = 1); u <- u/normVec(u)
  v <- stats::runif(dimX[2], min = -1, max = 1); v <- v/normVec(v)
  w <- stats::runif(dimX[3], min = -1, max = 1); w <- w/normVec(w)
  
  # initialize smoothing parameters: start with minimal values (almost no smoothing)
  alphaV <- min(alphaRange$v)
  alphaW <- min(alphaRange$w)
  
  # eigendecomposition of penalty matrices
  eV <- eigen(penMat$v)
  eW <- eigen(penMat$w)
  
  # digonal matrices
  Iv <- diag(dimX[2])
  Iw <- diag(dimX[3])
  
  # initialize results
  d <- rep(NA, K)
  
  U <- matrix(NA, ncol = K, nrow = dimX[1])
  V <- matrix(NA, ncol = K, nrow = dimX[2])
  W <- matrix(NA, ncol = K, nrow = dimX[3])
  
  for(k in seq_len(K))
  {
    if(verbose)
      cat("\nk = ", k, "\n")
    
    # initialize "old" versions
    uOld <- 0*u
    vOld <- 0*v
    wOld <- 0*w
    tolOld <- tol
    
    # number of iterations
    iter <- 0
    
    # repeat until convergence
    while(any(c(normVec(u - uOld)/normVec(u), normVec(v - vOld)/normVec(v), normVec(w - wOld)/normVec(w)) > tol))
    {
      # update u
      uOld <- u
      u <- ttv(X, list(v, w), dim = 2:3) /
        as.numeric(crossprod(v, v + alphaV*penMat$v %*%v )*crossprod(w, w + alphaW * penMat$w %*% w))
      
      # update v
      vOld <- v
      v <- solve(Iv + alphaV * penMat$v, ttv(X, list(u, w), dim = c(1,3)))/
        as.numeric(crossprod(u, u)*crossprod(w, w + alphaW * penMat$w %*% w))
      
      # update alphaV
      alphaV <- findAlphaVopt(alphaRange = alphaRange$v, data = X, u = u, w = w, 
                              alphaW = alphaW, OmegaW = penMat$w, 
                              GammaV = eV$vectors, lambdaV = eV$values)
         
      # update w
      wOld <- w
      w <- solve(Iw + alphaW * penMat$w, ttv(X, list(u, v), dim = 1:2)) /
        as.numeric(crossprod(u,  u)* crossprod(v, v + alphaV*penMat$v %*% v ))
      
      # update alphaW
      alphaW <- findAlphaWopt(alphaRange = alphaRange$w, data = X, u = u, v = v, 
                              alphaV = alphaV, OmegaV = penMat$v, 
                              GammaW = eW$vectors, lambdaW = eW$values)
      
      iter <- iter + 1
      
      if(iter > maxIter)
      {
        if((iter < 2*maxIter) & adaptTol) # try another maxIter with increased tolerance
        {
          tol <- 10*tol
        } 
        else 
        {
          uOld <- u
          vOld <- v
          wOld <- w
          warning("FCP-TPA algorithm did not converge; iteration ", k , " stopped.")
        } 
      }
    }
    
    if(verbose)
      cat("Absolute error:\n
          u: ", normVec(u - uOld), ", v: ", normVec(v - vOld), ", w: ",  normVec(w - wOld), ", alphaV: ", alphaV, ", alphaW: ", alphaW, "\n")
    
    # reset tolerance if necessary
    if(adaptTol & (iter >= maxIter) ) 
      tol <- tolOld
    
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

#' Check input of FCP TPA function
#' 
#' @return No return value, called for side effects
#' 
#' @keywords internal
check_FCP_TPA_input <- function(X, K, penMat, alphaRange, verbose, tol, maxIter, adaptTol)
{
  # X
  if(!is.array(X))
    stop("Parameter 'X' must be passed as an array.")
  else
  {
    if(length(dim(X)) != 3)
      stop("Parameter 'X' must have three dimensions.")
  }
  
  dimX <- dim(X)
  
  # K
  if(! (is.numeric(K) & length(K) == 1))
    stop("Parameter 'K' must be passed as a number.")
  
  # penMat
  if(!is.list(penMat))
    stop("Parameter 'penMat' must be passed as a list.")
  else
  {
    if(!identical(names(penMat), c("v", "w")))
      stop("Function FCP_TPA: penMat must be a list of matrices with entries 'v' and 'w'.")
    
    # penMat$v: correct dimensions and symmetric
    if(any(dim(penMat$v) != rep(dimX[2],2)))
      stop("Function FCP_TPA: the penalization matrix for dimension v must be ", dimX[2], " x ", dimX[2], ".")
    
    if(!isSymmetric(penMat$v))
      stop("Function FCP_TPA: the penalization matrix for dimension v must be symmetric.")
    
    # penMat$w: correct dimensions and symmetric
    if(any(dim(penMat$w) != rep(dimX[3],2)))
      stop("Function FCP_TPA: the penalization matrix for dimension w must be ", dimX[3], " x ", dimX[3], ".")
    
    if(!isSymmetric(penMat$w))
      stop("Function FCP_TPA: the penalization matrix for dimension w must be symmetric.")
  }
  
  # alphaRange
  if(!is.list(alphaRange))
    stop("Parameter 'alphaRange' must be passed as a list.")
  else
  {
    if(!identical(names(alphaRange), c("v", "w")))
      stop("Function FCP_TPA: alphaRange must be a list of vectors with entries 'v' and 'w'.")
    
    if(length(alphaRange$v) != 2)
      stop("Function FCP_TPA: alphaRange$v must be a vector of length 2.")
    
    if(any(alphaRange$v < 0))
      stop("Function FCP_TPA: Values for alphaV must not be negative.")
    
    if(length(alphaRange$w) != 2)
      stop("Function FCP_TPA: alphaRange$w must be a vector of length 2.")
    
    if(any(alphaRange$w < 0))
      stop("Function FCP_TPA: Values for alphaW must not be negative.")
  }
  
  # verbose
  if(!is.logical(verbose))
    stop("Parameter 'verbose' must be passed as a logical.")
  
  # tol
  if(! (is.numeric(tol) & length(tol) == 1))
    stop("Parameter 'tol' must be passed as a number.")
  
  # maxIter
  if(! (is.numeric(maxIter) & length(maxIter) == 1))
    stop("Parameter 'maxIter' must be passed as a number.")
  
  # adaptTol
  if(!is.logical(adaptTol))
    stop("Parameter 'adaptTol' must be passed as a logical.")
  
  return()
}

#' Calculate the euclidean norm of a vector
#' 
#' @param x The vector for which the norm is to be calculated
#' 
#' @return The euclidean norm of \code{x}.
#' 
#' @keywords internal
normVec <- function(x)
{
  sqrt(sum(x^2))
}

#' Find the optimal smoothing parameters in FCP_TPA using GCV
#' 
#' These functions find the optimal smoothing parameters \eqn{\alpha_v,
#' \alpha_w} for the two image directions (v and w) in the FCP_TPA algorithm
#' based on generalized cross-validation, which is nested in the tensor power
#' algorithm. Given a range of possible values of \eqn{\alpha_v} (or
#' \eqn{\alpha_w}, respectively), the optimum is found by optimizing the GCV
#' criterion using the function \code{\link[stats]{optimize}}.
#' 
#' @param alphaRange A numeric vector with two elements, containing the minimal 
#'   and maximal value for the smoothing parameter that is to be optimized.
#' @param data The tensor containing the data, an array of dimensions \code{N x 
#'   S1 x S2}.
#' @param u,v,w The current value of the eigenvectors \eqn{u_k, v_k, w_k} (not 
#'   normalized) of dimensions \code{N}, \code{S1} and \code{S2}.
#' @param alphaV,alphaW The current value of the smoothing parameter for the
#'   other image direction (\eqn{\alpha_w} for \code{findAlphaVopt} and
#'   \eqn{\alpha_v} for \code{findAlphaWopt}), which is kept as fixed.
#' @param OmegaV, OmegaW A matrix of dimension \code{S1 x S1} (\code{OmegaV} in  \code{findAlphaWopt}) or \code{S2 x S2} (\code{OmegaW} in  \code{findAlphaVopt}), the penalty matrix for 
#'  other image direction.
#' @param GammaV,GammaW A matrix of dimension \code{S1 x S1} (\code{GammaV} in  \code{findAlphaVopt}) or \code{S2 x S2} (\code{GammaW} in  \code{findAlphaWopt}), containing the 
#'   eigenvectors of the penalty matrix for the image direction for which the optimal smoothing parameter is to be found.
#' @param lambdaV, lambdaW A numeric vector of length  \code{S1}(\code{lambdaV} in  \code{findAlphaVopt}) or \code{S2} (\code{lambdaW} in  \code{findAlphaWopt}), containing the 
#'   eigenvalues of the penalty matrix for the image direction for which the optimal smoothing parameter is to be found.
#'   
#' @return  The optimal \eqn{\alpha_v} (or \eqn{\alpha_w}, respectively), found by optimizing the GCV criterion 
#'   within the given range of possible values.
#'   
#' @references G. I. Allen (2013), "Multi-way Functional Principal Components 
#'   Analysis", IEEE International Workshop on Computational Advances in 
#'   Multi-Sensor Adaptive Processing.
#' @references J. Z. Huang, H. Shen and A. Buja (2009), "The Analysis of Two-Way
#'   Functional Data Using Two-Way Regularized Singular Value Decomposition". 
#'   Journal of the American Statistical Association, Vol. 104, No. 488, 1609 --
#'   1620.
#'   
#' @keywords internal
#'   
#' @seealso \code{\link{FCP_TPA}}, \code{\link{gcv}}
#' 
#' @importFrom stats optimize
findAlphaVopt <- function(alphaRange, data, u, w, alphaW, OmegaW, GammaV, lambdaV)
{
  z <-  crossprod(GammaV, ttv(data, list(u,w), c(1,3))) / (normVec(u) * normVec(w))
  wOw <- as.numeric(crossprod(w, OmegaW %*% w)) # without as.numeric, this is a 1x1 matrix
  eta <- 1/(1 + alphaW * wOw/normVec(w))
  
  res <- stats::optimize(f=gcv, interval=c(min(alphaRange), max(alphaRange)), 
                 n = length(lambdaV), z = z, eta = eta, lambda = lambdaV)$minimum
 
 return(res)
}

#' @describeIn findAlphaVopt
#' 
#' @importFrom stats optimize
findAlphaWopt <- function(alphaRange, data, u, v, alphaV, OmegaV, GammaW, lambdaW)
{
  z <-  crossprod(GammaW, ttv(data, list(u,v), c(1,2))) / (normVec(u) * normVec(v))
  vOv <- as.numeric(crossprod(v, OmegaV %*% v))
  eta <- 1/(1 + alphaV * vOv/normVec(v)) # without as.numeric, this is a 1x1 matrix

  res <- stats::optimize(f=gcv, interval=c(min(alphaRange), max(alphaRange)), 
                  n = length(lambdaW), z = z, eta = eta, lambda = lambdaW)$minimum
  
  return(res)
}

#' Generalized cross-validation for the FCP-TPA algorithm
#' 
#' These function calculates the generalized cross-validation criterion for the 
#' smoothing parameters \eqn{\alpha_v} or \eqn{\alpha_w} that are used in the 
#' \code{\link{FCP_TPA}} algorithm. As the criterion is symmetric in \eqn{v} and
#' \eqn{w}, this function implements a generic criterion, which is called by 
#' \code{\link{findAlphaVopt}}, \code{\link{findAlphaWopt}} with the correct values.
#' 
#' The criterion can be evaluated in a numerically efficient way, adopting the ideas in Huang, Shen and Buja (2008) to three-ways tensors. TODO!
#' 
#' @param alpha The current value of the smoothing parameter.
#' @param n The length of the dimension, for which the smoothing parameter is to
#'   be optimized.
#' @param z A vector of length \code{n}. See Details.
#' @param eta A vector of length \code{n}. See Details.
#' @param lambda A vector of length \code{n}, containing the eigenvalues of the
#'   penalty matrix corresponding the the current image direction.
#'   
#' @return The value of the GCV criterion.
#'   
#' @references G. I. Allen, "Multi-way Functional Principal Components 
#'   Analysis", IEEE International Workshop on Computational Advances in 
#'   Multi-Sensor Adaptive Processing, 2013.
#'   
#' @keywords internal
#'   
#' @seealso \code{\link{FCP_TPA}}
gcv <- function(alpha, n, z, eta, lambda)
{
  eS <- 1/(1 + alpha*lambda) # eigenvalues of S
  
  res <- sum((z* (1 - eta* eS))^2)/n
  res <- res/(1 - eta/n * sum(eS))^2
  
  return(res)
}