#' UMPCA: Uncorrelated Multilinear Principle Component Analysis
#'
#' This function implements the uncorrelated multilinear principal component
#' analysis for tensors of dimension 2, 3 or 4. The code is basically the same
#' as in the  MATLAB toolbox UMPCA by Haiping Lu (Link:
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/35432-uncorrelated-multilinear-principal-component-analysis-umpca},
#' see also references).
#'
#' @section Warning: As this algorithm aims more at uncorrelated features than
#'   at an optimal reconstruction of the data, hence it might give poor results
#'   when used for the univariate decomposition of images in MFPCA.
#'
#' @param TX The input training data in tensorial representation, the last mode
#'   is the sample mode. For \code{N}th-order tensor data, \code{TX} is of
#'   \code{(N+1)}th-order with the \code{(N+1)}-mode to be the sample mode.
#'   E.g., 30x20x10x100 for 100 samples of size 30x20x10.
#' @param numP The dimension of the projected vector, denoted as \eqn{P} in the
#'   paper. It is the number of elementary multilinear projections (EMPs) in
#'   tensor-to-vector projection.
#'
#' @return \item{Us}{The multilinear projection, consisting of \code{numP}
#'   (\eqn{P} in the paper) elementary multilinear projections (EMPs), each EMP
#'   is consisted of \code{N} vectors, one in each mode.} \item{TXmean}{The mean
#'   of the input training samples \code{TX}.} \item{odrIdx}{The ordering index
#'   of projected features in decreasing variance.}
#'
#' @references Haiping Lu, K.N. Plataniotis, and A.N. Venetsanopoulos,
#'   "Uncorrelated Multilinear Principal Component Analysis for Unsupervised
#'   Multilinear Subspace Learning", IEEE Transactions on Neural Networks, Vol.
#'   20, No. 11, Page: 1820-1836, Nov. 2009.
#'
#' @export UMPCA
#'
#' @examples
#' set.seed(12345)
#'
#'  # define "true" components
#'  a <- sin(seq(-pi, pi, length.out = 100))
#'  b <- exp(seq(-0.5, 1, length.out = 150))
#'
#'  # simulate tensor data
#'  X <- a %o% b %o% rnorm(80, sd = 0.5)
#'
#'  # estimate one component
#'  UMPCAres <- UMPCA(X, numP = 1)
#'
#'  # plot the results and compare to true values
#'  plot(UMPCAres$Us[[1]][,1])
#'  points(a/sqrt(sum(a^2)), pch = 20) # eigenvectors are defined only up to a sign change!
#'  legend("topright", legend = c("True", "Estimated"), pch = c(20, 1))
#'
#'  plot(UMPCAres$Us[[2]][,1])
#'  points(b/sqrt(sum(b^2)), pch = 20)
#'  legend("topleft", legend = c("True", "Estimated"), pch = c(20, 1))
UMPCA <- function(TX, numP)
{
  if(!is.array(TX))
    stop("Parameter 'TX' must be passed as an array.")
  if(! (is.numeric(numP) & length(numP) == 1))
    stop("Parameter 'numP' must be passed as a number.")
  
  # TX: (N+1)-dimensional tensor Tensor Sample Dimension x NumSamples
  N <- length(dim(TX)) - 1 # the order of samples
  IsTX <- dim(TX)  
  Is <- IsTX[seq_len(N)]  # the dimensions of the tensor
  
  #Please see Corollary 1 in the TNN paper: numP<= min(min(Is),M}
  if(numP > min(Is))
    numP <- min(Is)
  
  numSpl <- IsTX[N+1] # number of samples
  if(numSpl < numP)
    numP <- numSpl-1 # centering makes the samples to be dependent, hence rank lowered by 1
  
  #### Zero-Mean ####
  TXmean <- rowMeans(TX, dims = N) # the mean
  for(n in seq_len(numSpl))
    TX[,,n] <- TX[,,n] - TXmean # centering
  dim(TXmean) <- c(dim(TXmean), 1) # correct dimensions for output
  
  #### UMPCA parameters ####
  maxK <- 10 # maximum number of iterations, you can change this number
  Us <- vector("list", N) # N matrices
  for(n in seq_len(N))
    Us[[n]] <- matrix(NA, nrow = Is[n], ncol = numP)
  Us0 <- vector("list", N) # N vectors
  for(iP in seq_len(numP)) # get each EMP one by one
  {
    #Initialization
    for(n in seq_len(N))
    {
      if(iP == 1)
      {
        Un <- rep(1, Is[n])
        Un <- Un / sqrt(sum(Un^2))   
        Us0[[n]] <- Un  
      }
      Us[[n]][,iP] <- Us0[[n]]  
    }
    #End Initialization
    
    #Start iterations
    for(k in seq_len(maxK))
    {
      for(n in seq_len(N))
      {
        idx <- seq_len(N)[-n]
        Ypn <- ttv(TX, sapply(Us[idx], function(x){x[,iP]}, simplify = FALSE), idx)

        ST <- rep(0, Is[n])  
        for(i in seq_len(numSpl))
        {
          YDiff <- Ypn[,i]
          ST <- ST + tcrossprod(YDiff, YDiff)   #Within-class Scatter
        }
        
        #### #### ####
        if(iP > 1)
        {
          GYYG <- crossprod(Ypn %*% Gps) #equation (16) in the paper
          ThetaST <- (diag(1, nrow = Is[n]) - Ypn %*% Gps %*% solve(GYYG) %*% t(Gps) %*% t(Ypn)) #equation (15) in the paper
          ST <- ThetaST %*% ST  
          Un <- maxeig(ST)$x
        }
        else {
          Un <- maxeig(ST)$x
        }
        Un <- Un / sqrt(sum(Un^2))
        Us[[n]][,iP] <- Un  
      }
    }
    
    #### Projection ####
    gp <- ttv(TX, sapply(Us, function(x){x[,iP]}, simplify = FALSE), seq_len(N))
    if (iP == 1)
      Gps <- gp  
    else
      Gps <- cbind(Gps, gp)
    #### End Projection ####
  }
  
  #### Sort by Variance ####
  vecYps <- t(Gps)  #vecYps contains the feature vectors for training data
  Ymean <- rowMeans(vecYps)  #Should be zero
  TVars <- diag( tcrossprod(vecYps, vecYps) )  #Calculate variance
  odrIdx <- order(TVars, decreasing = TRUE)  
  odrIdx <- odrIdx[seq_len(numP)] # take the first numP
  
  return(list(Us = Us, TXmean = TXmean, odrIdx = odrIdx))
}


#' Compute the largest eigenvalue and associated eigenvector of a matrix A using the power method
#' 
#' MATLAB Function: Copyright 1999 by Todd K. Moon
#' 
#' @param A Matrix whose eigenvalue is sought
#' 
#' @return 
#' \item{lambda}{Largest eigenvalue}
#' \item{x}{corresponding eigenvector}
#' 
#' @importFrom stats runif
#' 
#' @keywords internal
maxeig <- function(A)
{
  d <- dim(A)
  if(d[1] != d[2])
    stop("Input matrix must be quadratic!")
  n <- d[1]
  x <- rep(0, n)
   x[1] <- 1 
  # make sure x is not an eigenvector of A (for a possibly small eigenvalue)
  while( isTRUE( all.equal(as.vector(A %*% x), as.numeric(t(x) %*% A %*% x) * x) ))
  {
    x <- stats::runif(n)
    x <- x/normVec(n)
  }             
  lambda <- 1
  lambdaold <- 0
  maxItr <- 300
  iItr <- 1
  while(abs(lambda - lambdaold) > .Machine$double.eps & iItr < maxItr)
  {
    lambdaold <- lambda 
    z <-  A %*% x
    x <-  z / sqrt(sum(z^2))
    lambda <- t(x) %*% A %*% x
    iItr <- iItr+1 
  }  
  
  return(list(lambda = lambda, x = x)) 
}

#' Internal function for the Tensor times Vector calculation
#' 
#' @seealso \code{\link{ttv}}
#'  
#' @importFrom plyr alply
#' 
#' @keywords internal
ttvCalculation <- function(A, v, dim)
{
  res <- Reduce('+', mapply(function(arr,vec){arr*vec}, plyr::alply(A, dim), unlist(v), SIMPLIFY = FALSE))
  
  return(res)
}

#' Tensor times vector calculation
#' 
#' Functionality adapted from the MATLAB tensor toolbox 
#' (\url{http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html}).
#' 
#' Let \code{A} be a tensor with dimensions \eqn{d_1 \times d_2 \times \ldots 
#' \times d_p}{d_1 x d_2 x \ldots x d_p} and let \code{v} be a vector of length 
#' \eqn{d_i}. Then the tensor-vector-product along the \eqn{i}-th dimension is 
#' defined as \deqn{B_{j_1 \ldots j_{i-1}j_{i+1} \ldots j_d} = \sum_{i=1}^{d_i} 
#' A_{j_1 \ldots j_{i-1} i j_{i+1} \ldots j_d} \cdot v_i.}{ B[j_1, \ldots 
#' ,j_{i-1},j_{i+1},\ldots,j_d] = \sum A[j_1, \ldots, j_{i-1}, i, j_{i+1}, 
#' \ldots, j_d]  v[i].} It can hence be seen as a generalization of the 
#' matrix-vector product.
#' 
#' The tensor-vector-product along several dimensions between a tensor \code{A} 
#' and multiple vectors \code{v_1,\ldots,v_k} (\eqn{k \le p}) is defined as a 
#' series of consecutive tensor-vector-product along the different dimensions.
#' For consistency, the multiplications are calculated from the dimension of the
#' highest order to the lowest.
#' 
#' @param A An array.
#' @param v A list of the same length as \code{dim}.
#' @param dim A vector specifying the dimensions for the multiplication.
#'   
#' @return An array, the result of the multiplication.
#'   
#' @seealso \code{\link{UMPCA}}
#'   
#' @references B. W. Bader and T. G. Kolda. Algorithm 862: MATLAB tensor classes
#'   for fast algorithm prototyping, ACM Transactions on Mathematical Software 
#'   32(4):635-653, December 2006.
#'   
#' @export ttv
#' 
#' @examples
#' # create a three-mode tensor
#' a1 <- seq(0,1, length.out = 10)
#' a2 <- seq(-1,1, length.out = 20)
#' a3 <- seq(-pi, pi, length.out = 15)
#' A <-a1 %o% a2 %o% a3
#' dim(A)
#' 
#' # multiply along different dimensions
#' dim(ttv(A = A, v = list(rnorm(10)), dim = 1))
#' dim(ttv(A = A, v = list(rnorm(20)), dim = 2))
#' dim(ttv(A = A, v = list(rnorm(15)), dim = 3))
#' 
#' # multiply along more than one dimension
#' length(ttv(A = A, v = list(rnorm(10), rnorm(15)), dim = c(1,3)))
ttv <- function(A, v, dim)
{
  if(!is.array(A))
    stop("Parameter 'A' must be an array.")
  if(!is.list(v))
    stop("Parameter 'v' must be passed as a list.")
  if(!is.numeric(dim))
    stop("Parameter 'dim' must be passed as a vector of numerics.")
  
  # check input arguments
  if(any(dim(A)[dim] != vapply(v, FUN = length, FUN.VALUE = 0)))
    stop("A and v have wrong dimensions!")
  
  if(length(dim) != length(v))
    stop("The parameters 'dim' and 'v' must have the same length!")
  
  # reorder dimensions if necessary
  if(length(dim) > 1)
  {
    dimOrd <- order(dim, decreasing = TRUE) # order dimensions in descending order
    dim <- dim[dimOrd] # reorder dim argument
    v <- v[dimOrd] # reorder v list
  }
  
  for(d in seq_len(length(dim)))
    A <- ttvCalculation(A, v[d], dim[d])
    
  return(A)
}
