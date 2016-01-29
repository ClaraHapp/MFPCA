#' The functional CP-TPA algorithm
#' 
#' @param X The data tensor of dimensions \code{N x S1 x S2}
#' @param K The number of eigentensors to be calculated
#' @param penMat A list of three roughness penalization matrices, one for each direction (u, v & w).
#' @param alphaGrid A list of three grids for smoothing parameters along each direction. For the observations (first diection), only the minimum is used
#'   
#' @references G. Allen Multi-Way Functional Principal Component Analysis
FCP_TPA <- function(X, K, penMat, alphaGrid)
{
  dimX <- dim(X)
  
  # initialize the norm-one vectors u,v,w randomly
  u <- runif(dimX[1], min = -1, max = 1); u <- u/normVec(u)
  v <- runif(dimX[2], min = -1, max = 1); v <- v/normVec(v)
  w <- runif(dimX[3], min = -1, max = 1); w <- w/normVec(w)
  
  # initialize smoothness parameters and smoothing matrices
  alphaU <- min(alphaGrid$u); if(length(alphaGrid$u) > 1) warning("Set smoothing parameter alphaU (smoothing over observations) to minimum of given values. No optimization.") # not updated
  alphaV <- min(alphaGrid$v) # start with little smoothing
  alphaW <- min(alphaGrid$w)
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
    # cat("\nk = ", k, "\n")
    
    # initialize "old" versions
    uOld <- 0*u
    vOld <- 0*v
    wOld <- 0*w
    
    # repeat until convergence
    while(any(c(normVec(u - uOld)/normVec(u), normVec(v - vOld)/normVec(v), normVec(w - wOld)/normVec(w)) > sqrt(.Machine$double.eps)))
    {
      
      # update u (assume alpha_u = min(alphagrid$u) fixed, almost no smoothing over observations)
      uOld <- u
      u <- solve(Iu + alphaU * penMat$u, ttv(X, list(v, w), dim = 2:3)) /
        as.numeric(crossprod(v, v + alphaV*penMat$v %*%v )*crossprod(w, w + alphaW * penMat$w %*% w))#as.numeric(crossprod(v, solve(Sv, v))* crossprod(w, solve(Sw, w)))
      
      # update v
      vOld <- v
      v <- solve(Iv + alphaV * penMat$v, ttv(X, list(u, w), dim = c(1,3)))/
        as.numeric(crossprod(u, u + alphaU * penMat$u %*% u)*crossprod(w, w + alphaW * penMat$w %*% w))
      
      # optimize alphaV according to gcv on a given grid of values
      alphaV <- alphaGrid$v[which.min(sapply(X = alphaGrid$v, FUN = gcvV, data = X, u = u, v = v, w = w, penMat = penMat, alphaW = alphaW))]
      
      # update w
      wOld <- w
      w <- solve(Iw + alphaW * penMat$w, ttv(X, list(u, v), dim = 1:2)) /
        as.numeric(crossprod(u,  u + alphaU * penMat$u %*% u)* crossprod(v, v + alphaV*penMat$v %*% v ))
      
      # update alphaW according to gcv on a given grid of values
      alphaW <- alphaGrid$w[which.min(sapply(X = alphaGrid$w, FUN = gcvW, data = X, u = u, v=  v, w = w, penMat = penMat, alphaV = alphaV))]
      
      #cat("u: ", normVec(u - uOld), ", v: ", normVec(v - vOld), "w: ",  normVec(w - wOld), "alphaV: ", alphaV, "alphaW: ", alphaW, "\n")
    }
    
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


normVec <- function(x) sqrt(sum(x^2))

gcvV <- function(alphaV, data, u, v, w, penMat, alphaW)
{
  m <- length(v)
  
  res <- 1/m* normVec(ttv(data, list(u,w), c(1,3))/(normVec(u) * normVec(w))^2 - v)^2 / 
    (1 - normVec(w)^2/m * sum(1/eigen(diag(m) + alphaV*penMat$v)$values) /as.numeric(crossprod(w, w + alphaW * penMat$w %*% w)))^2  
  
  return(res)
}

gcvW <- function(alphaW, data, u, v, w, penMat, alphaV)
{
  p <- length(w)
  
  res <- 1/p* normVec(ttv(data, list(u,v), c(1,2))/(normVec(u) * normVec(v))^2 - w)^2 / 
    (1 - normVec(v)^2/p * sum(1 / eigen(diag(p) + alphaW * penMat$w)$values ) /as.numeric(crossprod(v, v + alphaV*penMat$v %*%v )))^2  
  
  return(res)
}
