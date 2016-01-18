#' UMPCA: Uncorrelated Multilinear Principle Component Analysis
#' 
#' This function implements the uncorrelated multilinear principal component 
#' analysis for tensors of dimension 2, 3 or 4. The code is basically the same as
#' in the  MATLAB toolbox
#' \url{UMPCA}{http://www.mathworks.com/matlabcentral/fileexchange/35432-uncorrelated-multilinear-principal-component-analysis--umpca-}
#' by Haiping Lu (see references).
#' 
#' @param TX The input training data in tensorial representation, the last mode
#'        is the sample mode. For Nth-order tensor data, \code{TX} is of 
#'        (N+1)th-order with the (N+1)-mode to be the sample mode.
#'        E.g., 30x20x10x100 for 100 samples of size 30x20x10#'
#' @param numP The dimension of the projected vector, denoted as \eqn{P} in the
#'          paper. It is the number of elementary multilinear projections 
#'          (EMPs) in tensor-to-vector projection.
#'          
#'  @return \item{}
#'    \item{Us}{The multilinear projection, consiting of \code{numP} (\eqn{P} in the paper) 
#'        elementary multilinear projections (EMPs), each EMP is consisted
#'        of \eqn{N} vectors, one in each mode.}
#'  \item{TXmean}{The mean of the input training samples \code{TX}.}
#'  \item{odrIdx}{The ordering index of projected features in decreasing  
#'            variance.}
#' 
#' @references 
#'  Haiping Lu, K.N. Plataniotis, and A.N. Venetsanopoulos,
#'   "Uncorrelated Multilinear Principal Component Analysis for Unsupervised Multilinear Subspace Learning",
#'    IEEE Transactions on Neural Networks,
#'   Vol. 20, No. 11, Page: 1820-1836, Nov. 2009.
UMPCA <- function(TX, numP)
{
  # TX: (N+1)-dimensional tensor Tensor Sample Dimension x NumSamples
  N <- length(dim(TX)) - 1 # the order of samples
  IsTX <- dim(TX)  
  Is <- IsTX[1:N]  # the dimensions of the tensor
  
  #Please see Corollary 1 in the TNN paper: numP<= min(min(Is),M}
  if(numP > min(Is))
    numP <- min(Is)
  
  numSpl <- IsTX[N+1] # number of samples
  if(numSpl < numP)
    numP <- numSpl-1 # centering makes the samples to be dependent, hence rank lowered by 1
  
  #### Zero-Mean ####
  TXmean <- rowMeans(TX, dims = N) # the mean
  for(n in 1:numSpl)
    TX[,,n] <- TX[,,n] - TXmean # centering
  dim(TXmean) <- c(dim(TXmean), 1) # correct dimensions for output
  
  #### UMPCA parameters ####
  maxK <- 10 # maximum number of iterations, you can change this number
  Us <- vector("list", N) # N matrices
  for(n in 1:N)
    Us[[n]] <- matrix(NA, nrow = Is[n], ncol = numP)
  Us0 <- vector("list", N) # N vectors
  for(iP in 1:numP) # get each EMP one by one
  {
    #Initialization
    for(n in 1:N)
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
    for(k in 1:maxK)
    {
      for(n in 1:N)
      {
        idx <- (1:N)[-n]
        Ypn <- ttv(TX, sapply(Us[idx], function(x){x[,iP]}, simplify = FALSE), idx)

        ST <- rep(0, Is[n])  
        for(i in 1:numSpl)
        {
          YDiff <- Ypn[,i]
          ST <- ST + YDiff %*% t(YDiff)   #Within-class Scatter
        }
        
        #### #### ####
        if(iP > 1)
        {
          GYYG <- t(Gps) %*% t(Ypn) %*% Ypn %*% Gps #equation (16) in the paper
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
    gp <- ttv(TX, sapply(Us, function(x){x[,iP]}, simplify = FALSE), 1:N)  
    if (iP == 1)
      Gps <- gp  
    else
      Gps <- cbind(Gps, gp)
    #### End Projection ####
  }
  
  #### Sort by Variance ####
  vecYps <- t(Gps)  #vecYps contains the feature vectors for training data
  Ymean <- rowMeans(vecYps)  #Should be zero
  TVars <- diag(vecYps %*% t(vecYps))  #Calculate variance
  odrIdx <- order(TVars, decreasing = TRUE)  
  odrIdx <- odrIdx[1:numP] # take the first numP
  
  return(list(Us = Us, TXmean = TXmean, odrIdx = odrIdx))
}


#' Compute the largest eigenvalue and associated eigenvector of a matrix A using the power method
#' 
#' Matlab Function: Copyright 1999 by Todd K. Moon
#' 
#' @param A Matrix whose eigenvalue is sought
#' 
#' @return 
#' \item{lambda}{Largest eigenvalue}
#' \item{x}{corresponding eigenvector}
#' 
#' @keywords internal
maxeig <- function(A)
{
  d <- dim(A)
  if(d[1] != d[2])
    stop("Input matrix must be quadratic!")
  n <- d[1]
  x <- rep(0, n)
  x[1] <- 1  # assumed to be not orthogonal to the first eigenvector
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
#' @seealso ttv
#' 
#' @keywords internal
ttvCalculation <- function(A, v, dim)
{
  res <- Reduce('+', mapply(function(arr,vec){arr*vec}, plyr::alply(A, dim), unlist(v), SIMPLIFY = FALSE))
  
  return(res)
}

#' Tensor times vector calculation
#' 
#' Functionality adapted from MATLAB tensor toolbox 
#' 
#' @param A An array
#' @param v A list of the same length as dim
#' @param dim A vector specifying the dimensions for the multiplication
#' 
#' @return An array, the result of the multiplication
ttv <- function(A, v, dim)
{
  # check input arguments
  if(any(dim(A)[dim] != sapply(v, length)))
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
  
  for(d in 1:length(dim))
    A <- ttvCalculation(A, v[d], dim[d])
    
  return(A)
}
