#' Univariate basis decomposition
#'
#' This function calculates a univariate basis decomposition for a (univariate)
#' functional data object.
#'
#' Functional data \eqn{X_i(t)} can often be approximated by a linear
#' combination of basis functions \eqn{b_k(t)} \deqn{X_i(t) = \sum_{k = 1}^K
#' \theta_{ik} b_k(t), i = 1, \ldots, N.} The basis functions may be
#' prespecified (such as spline basis functions or Fourier bases) or can be
#' estimated from the data (e.g. by functional principal component analysis) and
#' are the same for all observations \eqn{X_1(t), \ldots, X_n(t)}. The
#' coefficients (or scores) \eqn{\theta_{ik}} reflect the weight of each basis
#' function \eqn{b_k(t)} for the observed function \eqn{X_i(t)} and can be used
#' to characterize the individual observations.
#'
#' @param type A character string, specifying the basis for which the
#'   decomposition is to be calculated.
#' @param funDataObject A \code{funData} object, representing the (univariate) functional
#'   data samples.
#' @param ... Further parameters, passed to the function for the particular basis to
#'   use.
#'
#' @return \item{scores}{A matrix of scores (coefficients) for each observation
#'   based on the prespecified basis functions.} \item{B}{A matrix containing
#'   the scalar products of the basis functions. Can be \code{NULL} if the basis
#'   functions are orthonormal.} \item{ortho}{Logical. If \code{TRUE}, the basis
#'   functions are all orthonormal.} \item{functions}{A functional data object,
#'   representing the basis functions. Can be \code{NULL} if the basis functions
#'   are not estimated from the data, but have a predefined form. See Details.}
#'
#' @seealso \code{\link{MFPCA}}, \code{\link{univExpansion}}, \code{\link{fpcaBasis}}, \code{\link{splineBasis1D}},
#'   \code{\link{splineBasis1Dpen}}, \code{\link{splineBasis2D}}, \code{\link{splineBasis2Dpen}},
#'   \code{\link{umpcaBasis}}, \code{\link{fcptpaBasis}},
#'   \code{\link{dctBasis2D}}, \code{\link{dctBasis3D}}
#'   
#' @export univDecomp
#'   
#' @examples
#' # generate some data
#' dat <- simFunData(argvals = seq(0,1,0.01), M = 5, 
#'                   eFunType = "Poly", eValType = "linear", N = 100)$simData
#' 
#' # decompose the data in univariate functional principal components...
#' decFPCA <- univDecomp(type = "uFPCA", funDataObject = dat, npc = 5)
#' str(decFPCA)
#' 
#' # or in splines (penalized)
#' decSplines <- univDecomp(type = "splines1Dpen", funDataObject = dat) # use mgcv's default params
#' str(decSplines)
univDecomp <- function(type, funDataObject, ...)
{
  params <- as.list(match.call()) # get all arguments
  params$funDataObject <- funDataObject # add funDataObject (-> make sure is evaluated in correct env.)
  
  # check if type and data are of correct type
  if(is.null(params$type))
    stop("univDecomp: must specify 'type'.")
  
  if(!inherits(params$type, "character"))
    stop("univDecomp: 'type' must be of class character.")
  
  if(is.null(params$funDataObject))
    stop("univDecomp: must specify 'funDataObject'.")
  
  if(class(params$funDataObject) != "funData")
    stop("univDecomp: 'funDataObject' must be of class funData.")
  
  # delete function call and type information in params
  params[[1]] <- NULL
  params$type <- NULL  
  
  res <- switch(type,
                "uFPCA" = do.call(fpcaBasis, params),
                "UMPCA" = do.call(umpcaBasis, params),
                "FCP_TPA" = do.call(fcptpaBasis, params),
                "splines1D" = do.call(splineBasis1D, params),
                "splines1Dpen" = do.call(splineBasis1Dpen, params),
                "splines2D" = do.call(splineBasis2D, params),
                "splines2Dpen" = do.call(splineBasis2Dpen, params),
                "DCT2D" = do.call(dctBasis2D, params),
                "DCT3D" = do.call(dctBasis3D, params),
                stop("Univariate Decomposition for 'type' = ", type, " not defined!")
  )
  
  if(res$ortho == FALSE & is.null(res$B))
    stop("UnivDecomp: must provide integral matrix B for non-orthonormal basis functions.")
  
  return(res)
}

#' Calculate a functional principal component basis representation for 
#' functional data on one-dimensional domains
#' 
#' This function calculates a functional principal component basis
#' representation for functional data on one-dimensional domains. The FPCA is
#' calculated via the \code{\link{PACE}} function, which is built on
#' \link[refund]{fpca.sc} in the \pkg{refund} package.
#' 
#' @param funDataObject An object of class \code{\link[funData]{funData}} 
#'   containing the observed functional data samples and for which the FPCA is 
#'   to be calculated.
#' @param nbasis An integer, representing the number of  B-spline basis 
#'   functions used for estimation of the mean function and bivariate smoothing 
#'   of the covariance surface. Defaults to \code{10} (cf. 
#'   \code{\link[refund]{fpca.sc}}).
#' @param pve A numeric value between 0 and 1, the proportion of variance 
#'   explained: used to choose the number of principal components. Defaults to 
#'   \code{0.99} (cf. \code{\link[refund]{fpca.sc}}).
#' @param npc An integer, giving a prespecified value for the number of 
#'   principal components. Defaults to \code{NULL}. If given, this overrides 
#'   \code{pve} (cf. \code{\link[refund]{fpca.sc}}).
#' @param makePD Logical: should positive definiteness be enforced for the 
#'   covariance surface estimate? Defaults to \code{FALSE} (cf. 
#'   \code{\link[refund]{fpca.sc}}).
#'   
#' @return \item{scores}{A matrix of scores (coefficients) with dimension 
#'   \code{N x K}, reflecting the weights for each principal component in each 
#'   observation, where \code{N} is the number of observations in
#'   \code{funDataObject} and \code{K} is the number of functional principal
#'   components.} \item{ortho}{Logical, set to \code{TRUE}, as basis functions
#'   are orthonormal.} \item{functions}{A functional data object, representing
#'   the functional principal component basis functions.} \item{meanFunction}{The smoothed mean function.}
#'   
#' @seealso \code{\link{univDecomp}}, \code{\link{PACE}}
#' 
#' @export fpcaBasis
#' 
#' @examples
#' # simulate N = 100 observations of functional data based on polynomial eigenfunctions on [0,1]
#' sim <- simFunData(argvals = seq(0,1,0.01), M = 5, eFunType = "Poly", eValType = "linear", N = 100)
#' 
#' # estimate the first 5 functional principal components from the data
#' fpca <- fpcaBasis(sim$simData, npc = 5)
#' 
#' oldpar <- par(no.readonly = TRUE)
#' par(mfrow = c(1,2))
#' plot(sim$trueFuns, obs = 1:5, main = "True eigenfunctions")
#' plot(fpca$functions, main = "Estimated eigenfunctions")
#' 
#' # Flip if necessary
#' plot(sim$trueFuns, obs = 1:5, main = "True eigenfunctions")
#' plot(flipFuns(extractObs(sim$trueFuns, obs = 1:5), fpca$functions),
#'      main = "Estimated eigenfunctions\n(flipped)")
#' 
#' par(oldpar)
fpcaBasis <- function(funDataObject, nbasis = 10, pve = 0.99, npc = NULL, makePD = FALSE)
{
  FPCA <- PACE(funDataObject, predData = NULL, nbasis, pve, npc, makePD)
  
  return(list(scores = FPCA$scores,
              ortho = TRUE,
              functions = FPCA$functions,
              meanFunction = FPCA$mu
  ))
}

#' Calculate an uncorrelated multilinear principal component basis 
#' representation for functional data on two-dimensional domains
#' 
#' This function calculates an uncorrelated multilinear principal component 
#' analyis (UMPCA) representation for functional data on two-dimensional 
#' domains. In this case, the data can be interpreted as images with \code{S1 x 
#' S2} pixels (assuming \code{nObsPoints(funDataObject) = (S1, S2)}), i.e. the 
#' total observed data are represented as third order tensor of dimension 
#' \code{N x S1 x S2}.  The UMPCA of a tensor of this kind is calculated via the
#' \link{UMPCA} function, which is an \code{R}-version of the analogous 
#' functions in the \code{UMPCA} MATLAB toolbox by Haiping Lu (Link: 
#' \url{http://www.mathworks.com/matlabcentral/fileexchange/35432}, see also 
#' references).
#' 
#' @section Warning: As this algorithm aims more at uncorrelated features than 
#'   at an optimal reconstruction of the data, hence it might give poor results 
#'   when used for the univariate decomposition of images in MFPCA. The function
#'   therefore throws a warning.
#'   
#' @param funDataObject An object of class \code{\link[funData]{funData}} 
#'   containing the observed functional data samples (here: images) for which 
#'   the UMPCA is to be calculated.
#' @param npc An integer, giving the number of principal components to be 
#'   calculated.
#'   
#' @return \item{scores}{A matrix of scores (coefficients) with dimension 
#'   \code{N x k}, reflecting the weight of each principal component in each 
#'   observation.}  \item{B}{A matrix containing the scalar product of all pairs
#'   of basis functions.} \item{ortho}{Logical, set to \code{FALSE}, as basis 
#'   functions are not orthonormal.} \item{functions}{A functional data object, 
#'   representing the functional principal component basis functions.}
#'   
#' @seealso \code{\link{univDecomp}}
#'   
#' @references Haiping Lu, K.N. Plataniotis, and A.N. Venetsanopoulos, 
#'   "Uncorrelated Multilinear Principal Component Analysis for Unsupervised 
#'   Multilinear Subspace Learning", IEEE Transactions on Neural Networks, Vol. 
#'   20, No. 11, Page: 1820-1836, Nov. 2009.
#'   
#' @export umpcaBasis
#' 
#' @examples
#' # simulate image data for N = 100 observations
#' N <- 100
#' b1 <- eFun(seq(0,1,0.01), M = 7, type = "Poly")
#' b2 <- eFun(seq(-pi, pi, 0.03), M = 8, type = "Fourier")
#' b <- tensorProduct(b1,b2) # 2D basis functions
#' scores <- matrix(rnorm(N*56), nrow = N)
#' 
#' # calculate observations (= linear combination of basis functions)
#' f <- expandBasisFunction(scores = scores, functions = b)
#' 
#' # calculate basis functions based on UMPCA algorithm (needs some time)
#' \donttest{
#' # throws warning as the function aims more at  uncorrelated features than at
#' # optimal data reconstruction (see help)
#' umpca <- umpcaBasis(f, npc = 5) 
#' 
#' oldpar <- par(no.readonly = TRUE)
#' 
#' for(i in 1:5) # plot all 5 basis functions
#' plot(umpca$functions, obs = i, main = paste("Basis function", i)) # plot first basis function
#' 
#' par(oldpar)}
umpcaBasis <- function(funDataObject, npc)
{
  if(dimSupp(funDataObject) != 2)
    stop("UMPCA is implemented for (2D) image data only!")
  
  # throw warning
  warning("The UMPCA algorithm aims more at uncorrelated features than 
    at an optimal reconstruction of the data. It might hence give poor results 
    when used for the univariate decomposition of images in MFPCA.")
  
  # calculate UMPCA
  # permute observed data s.t. the observations are saved in the last dimension
  UMPCAres <- UMPCA(aperm(funDataObject@X, c(1+1:dimSupp(funDataObject), 1)), numP = npc)
  
  # calculcate eigenfunctions  
  eigenImages <- array(NA, c(npc, nObsPoints(funDataObject)))
  for(i in 1:npc)
    eigenImages[UMPCAres$odrIdx[i],,] <- tcrossprod(UMPCAres$Us[[1]][,i], UMPCAres$Us[[2]][,i])
  
  eigenFunctions <- funData(argvals = funDataObject@argvals, X = eigenImages)
  
  # calculate scores
  obsCent <- sweep(aperm(funDataObject@X, c(1+1:dimSupp(funDataObject), 1)), 1:2, UMPCAres$TXmean[,,1]) # demean images
  scores <- array(NA, c(nObs(funDataObject), npc))
  
  for(i in 1:npc)
    scores[,UMPCAres$odrIdx[i]] <- ttv(obsCent, sapply(UMPCAres$Us, function(x){x[,i]}, simplify = FALSE), 1:2)
  
  return(list(scores = scores,
              B = calcBasisIntegrals(eigenImages, dimSupp(funDataObject), funDataObject@argvals),
              ortho = FALSE,
              functions = eigenFunctions))
}


# Recursion for difference operator matrix
# taken from fpca.ssvd function in refund
makeDiffOp <- function(degree, dim){
  if(degree==0){
    return(diag(dim))
  } else {
    return(diff(makeDiffOp(degree-1, dim)))
  }
}


#' Calculate a smooth PCA representaion for functional data on two-dimensional 
#' domains
#' 
#' This function calculates a smooth PCA representation based on the FCP_TPA 
#' algorithm (see References) for functional data on two-dimensional domains. In
#' this case, the data can be interpreted as images with \code{S1 x S2} pixels 
#' (assuming \code{nObsPoints(funDataObject) = (S1, S2)}), i.e. the total data 
#' for \code{N} observations can be represented as third order tensor of 
#' dimension \code{N x S1 x S2}.
#' 
#' The smooth PCA of the tensor data is calculated via the \code{\link{FCP_TPA}}
#' function. Smoothness is induced by difference penalty matrices for both 
#' directions of the images, weighted by smoothing parameters \eqn{\alpha_v, 
#' \alpha_w}. The resulting eigenvectors can be interpreted in terms of 
#' eigenfunctions and individual scores for each observation. See 
#' \code{\link{FCP_TPA}} for details.
#' 
#' @param funDataObject An object of class \code{\link[funData]{funData}} 
#'   containing the observed functional data samples (here: images) for which 
#'   the smooth PCA is to be calculated.
#' @param npc An integer, giving the number of principal components to be 
#'   calculated.
#' @param smoothingDegree A numeric vector of length 2, specifying the degree of
#'   the difference penalties inducing smoothness in both directions of the 
#'   image. Defaults to \code{2} for each direction (2nd differences).
#' @param alphaRange A list of length 2 with entries \code{v} and \code{w} 
#'   containing the range of smoothness parameters to test for each direction.
#'   
#' @return \item{scores}{A matrix of scores (coefficients) with dimension 
#'   \code{N x K}, reflecting the weights for principal component in each 
#'   observation.}  \item{B}{A matrix containing the scalar product of all pairs
#'   of basis functions.} \item{ortho}{Logical, set to \code{FALSE}, as basis 
#'   functions are not orthonormal.} \item{functions}{A functional data object, 
#'   representing the functional principal component basis functions.}
#'   
#' @seealso \code{\link{univDecomp}}, \code{\link{FCP_TPA}}
#'   
#' @references G. I. Allen, "Multi-way Functional Principal Components 
#'   Analysis", In IEEE International Workshop on Computational Advances in 
#'   Multi-Sensor Adaptive Processing, 2013.
#'   
#' @export fcptpaBasis
#'   
#' @examples
#' # simulate image data for N = 100 observations
#' N <- 100
#' b1 <- eFun(seq(0,1,0.01), M = 7, type = "Poly")
#' b2 <- eFun(seq(-pi, pi, 0.03), M = 8, type = "Fourier")
#' b <- tensorProduct(b1,b2) # 2D basis functions
#' scores <- matrix(rnorm(N*56), nrow = N)
#' 
#' # calculate observations (= linear combination of basis functions)
#' f <- expandBasisFunction(scores = scores, functions = b)
#' 
#' # calculate basis functions based on FCP_TPA algorithm (needs some time)
#' \donttest{
#' fcptpa <- fcptpaBasis(f, npc = 5, alphaRange = list(v = c(1e-5, 1e5), w = c(1e-5, 1e5)))
#' 
#' oldpar <- par(no.readonly = TRUE)
#' 
#' for(i in 1:5) # plot all 5 basis functions
#' plot(fcptpa$functions, obs = i, main = paste("Basis function", i)) # plot first basis function
#' 
#' par(oldpar)}
fcptpaBasis <- function(funDataObject, npc, smoothingDegree = rep(2,2), alphaRange)
{
  if(dimSupp(funDataObject) != 2)
    stop("FCP_TPA is implemented for (2D) image data only!")
  
  d <- dim(funDataObject@X)
    
  # smoothing only along image directions
  Dv <- crossprod(makeDiffOp(degree = smoothingDegree[1], dim = d[2]))
  Dw <- crossprod(makeDiffOp(degree = smoothingDegree[2], dim = d[3]))
  
  pca <-  FCP_TPA(X = funDataObject@X, K = npc, penMat = list(v = Dv, w = Dw), alphaRange = alphaRange)
    
  # reconstruct eigenimages and scores from FCP_TPA result
  eigenImages <- array(NA, c(npc, d[-1]))
  
  for(i in 1:npc)
    eigenImages[i,,] <-  pca$V[,i] %o% pca$W[,i]
  
  scores <-  sweep(pca$U,MARGIN=2,pca$d,`*`)
  
  return(list(scores = scores,
              B = calcBasisIntegrals(eigenImages, dimSupp(funDataObject), funDataObject@argvals),
              ortho = FALSE,
              functions = funData(argvals = funDataObject@argvals, X = eigenImages)))
}


#' Calculate a spline basis decomposition for functional data on one-dimensional
#' domains
#' 
#' These functions calculate a penalized or unpenalized spline basis 
#' decomposition for functional data on one-dimensional domains based on the 
#' \link[mgcv]{gam} function in the \pkg{mgcv} package.
#' 
#' @param funDataObject An object of class \code{\link[funData]{funData}} 
#'   containing the observed functional data samples and for which the basis 
#'   decomposition is calculated.
#' @param bs A character string, specifying the type of basis functions to be 
#'   used. Defaults to \code{"ps"} (B-spline functions). Please refer to 
#'   \code{\link[mgcv]{smooth.terms}} for a list of possible basis functions.
#' @param m A numeric, the order of the spline basis. Defaults to \code{NA}, 
#'   i.e. the order is chosen automatically. See  \code{\link[mgcv]{s}} for 
#'   details.
#' @param k A numeric, the number of basis functions used. Defaults to 
#'   \code{-1}, i.e. the number of basis functions is chosen automatically. See 
#'   \code{\link[mgcv]{s}} for details.
#' @param parallel Logical (only for \code{splineBasis1Dpen}. If \code{TRUE}, 
#'   the coefficients for the basis functions are calculated in parallel. The 
#'   implementation is based on the \code{\link[foreach]{foreach}} function and 
#'   requires a parallel backend that must be registered before. See 
#'   \code{\link[foreach]{foreach}} for details.
#'   
#' @return \item{scores}{A matrix of scores (coefficients) with dimension 
#'   \code{N x K}, reflecting the weights for each of the \code{K} basis 
#'   functions and for each of the \code{N} observations.} \item{B}{A matrix 
#'   containing the scalar product of all pairs of basis functions.} 
#'   \item{ortho}{Logical, set to \code{FALSE}, as basis functions are not 
#'   orthonormal.} \item{functions}{\code{NULL}, as basis functions are known} 
#'   \item{settings}{A list with entries \code{bs}, \code{m} and \code{k}, 
#'   giving the actual parameters used for generating the spline basis 
#'   functions.}
#'   
#' @seealso \code{\link{univDecomp}}, \code{\link[mgcv]{gam}},
#'   \code{\link[foreach]{foreach}}
#'   
#' @importFrom stats lm
#' @importFrom mgcv gam s
#'   
#' @export splineBasis1D
#' 
#' @examples
#' # generate some data
#' dat <- simFunData(argvals = seq(0,1,0.01), M = 5, 
#'                   eFunType = "Poly", eValType = "linear", N = 100)$simData
#'                   
#'  # calculate spline basis decomposition
#'  dataDec <- splineBasis1D(dat) # use mgcv's default parameters
#'  str(dataDec)
#'  
#'  # add some noise to the data
#'  noisyDat <- addError(dat, sd = 0.5)
#'  
#'  # calculate spline basis decomposition with penalization to reduce noise
#'  noisyDataDec <- splineBasis1Dpen(dat) # use mgcv's default parameters
#'  str(noisyDataDec)
#'  
#'  # check if noise has been filtered out by penalization
#'  all.equal(noisyDataDec$scores, dataDec$scores, check.attributes = FALSE)
#'  # -> have almost the same coefficients
splineBasis1D <- function(funDataObject, bs = "ps", m = NA, k = -1)
{
  if(dimSupp(funDataObject) != 1)
    stop("splines1D is implemented for 1D functional data only.")
  
  N <- nObs(funDataObject)
  
  x <- funDataObject@argvals[[1]]
  
  # spline design matrix via gam
  g <- mgcv::gam(funDataObject@X[1, ] ~ s(x, bs = bs, m = m, k = k), fit = FALSE)
  desMat <- g$X
  k <- g$smooth[[1]]$bs.dim
  m <- g$smooth[[1]]$p.order
  
  # weights via lm -> no penalization
  scores <- t(apply(funDataObject@X, 1, function(f, dM){stats::lm(f ~ dM - 1)$coef}, dM = desMat)) # design matrix already includes intercept!
  
  return(list(scores = scores,
              B = calcBasisIntegrals(t(desMat), 1, funDataObject@argvals),
              ortho = FALSE,
              functions = NULL,
              settings = list(bs = bs, k = k, m = m)
  ))
}

#' @rdname splineBasis1D
#'
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#' @importFrom mgcv gam model.matrix.gam s
#' 
#' @export splineBasis1Dpen
splineBasis1Dpen <- function(funDataObject, bs = "ps", m = NA, k = -1, parallel = FALSE)
{
  if(dimSupp(funDataObject) != 1)
    stop("splines1Dpen is implemented for 1D functional data only.")
  
  N <- nObs(funDataObject)
  
  x <- funDataObject@argvals[[1]]
  
  if(parallel)
  {
    scores <- foreach::foreach(i = 1:(N-1), .combine = "rbind")%dopar%{
      g <- mgcv::gam(funDataObject@X[i, ] ~ s(x, bs = bs, m = m, k = k), method = "REML")
      g$coef
    }
  }
  else
  {
    scores <- foreach::foreach(i = 1:(N-1), .combine = "rbind")%do%{
      g <- mgcv::gam(funDataObject@X[i, ] ~ s(x, bs = bs, m = m, k = k), method = "REML")
      g$coef
    }
  }
  
  # last extra to extract model matrix
  g <- mgcv::gam(funDataObject@X[N, ] ~ s(x, bs = bs, m = m, k = k), method = "REML")
  k <- g$smooth[[1]]$bs.dim
  m <- g$smooth[[1]]$p.order
  
  scores <- rbind(scores, g$coef)
  
  return(list(scores = scores,
              B = calcBasisIntegrals(t(mgcv::model.matrix.gam(g)), 1, funDataObject@argvals),
              ortho = FALSE,
              functions = NULL,
              settings = list(bs = bs, k = k, m = m)
  ))
}


#' Calculate a spline basis representation for functional data on 
#' two-dimensional domains
#' 
#' These functions calculate a penalized or unpenalized tensor product spline 
#' basis representation for functional data on two-dimensional domains based on 
#' the \code{\link[mgcv]{gam}}/\code{\link[mgcv]{bam}} functions in the
#' \pkg{mgcv} package. See Details.
#' 
#' If the basis representation is calculated without penalization 
#' (\code{splineBasis2D}), the coefficients are computed using the 
#' \code{\link[mgcv]{gam}} function from the \pkg{mgcv} package. In the case of 
#' penalization (\code{splineBasis2Dpen}), the function \code{\link[mgcv]{bam}}
#' (for large GAMs) is used instead.
#' 
#' @param funDataObject An object of class \code{\link[funData]{funData}} 
#'   containing the observed functional data samples and for which the basis 
#'   representation is calculated.
#' @param bs A vector of character strings (or a single character string), 
#'   specifying the type of basis functions to be used. Defaults to \code{"ps"} 
#'   (P-spline functions). Please refer to \code{\link[mgcv]{te}} for a list of 
#'   possible basis functions.
#' @param m A numeric vector (or a single number), the order of the spline 
#'   basis. Defaults to \code{NA}, i.e. the order is chosen automatically.  See 
#'   \code{\link[mgcv]{s}} for details.
#' @param k An numeric vector (or a single number), the number of basis 
#'   functions used.  Defaults to \code{-1}, i.e. the number of basis functions 
#'   is chosen automatically.   See  \code{\link[mgcv]{s}} for details.
#' @param parallel Logical (only for function \code{splineBasis2Dpen}). If 
#'   \code{TRUE}, the coefficients for the basis functions are calculated in 
#'   parallel. The implementation is based on the \code{\link[foreach]{foreach}}
#'   function and requires a parallel backend that must be registered before. 
#'   See \code{\link[foreach]{foreach}} for details.
#'   
#' @return \item{scores}{A matrix of scores (coefficients) with dimension 
#'   \code{N x K}, reflecting the weights for each basis function in each 
#'   observation, where \code{K} is the total number of basis functions used.} 
#'   \item{B}{A matrix containing the scalar product of all pairs of basis 
#'   functions.} \item{ortho}{Logical, set to \code{FALSE}, as basis functions 
#'   are not orthonormal.} \item{functions}{\code{NULL}, as basis functions are 
#'   known.} \item{settings}{A list with entries \code{bs}, \code{m} and 
#'   \code{k}, giving the actual parameters used for generating the spline basis
#'   functions.}
#'   
#' @seealso \code{\link{univDecomp}}, \code{\link{splineBasis1D}}, 
#'   \code{\link[mgcv]{gam}}, \code{\link[mgcv]{bam}}, 
#'   \code{\link[foreach]{foreach}}
#'   
#' @importFrom stats lm
#' @importFrom mgcv gam te
#'   
#' @export splineBasis2D
#' 
#' @examples
#' # simulate image data for N = 100 observations
#' N <- 100
#' b1 <- eFun(seq(0,1,0.01), M = 7, type = "Poly")
#' b2 <- eFun(seq(-pi, pi, 0.03), M = 8, type = "Fourier")
#' b <- tensorProduct(b1,b2) # 2D basis functions
#' scores <- matrix(rnorm(N*56), nrow = N)
#' 
#' # calculate observations (= linear combination of basis functions)
#' dat <- expandBasisFunction(scores = scores, functions = b)
#' 
#' # calculate 2D spline basis decomposition (needs some time)
#' \donttest{
#' dataDec <- splineBasis2D(dat, k = c(5,5)) # use 5 basis functions in each direction
#' }
#' 
#' # add some noise to the data
#' noisyDat <- addError(dat, sd = 0.5)
#' 
#' # calculate 2D spline basis decomposition with penalization (needs A LOT more time)
#' \donttest{
#' noisyDataDec <- splineBasis2Dpen(noisyDat, k = c(5,5)) # use 5 basis functions in each direction
#' }
splineBasis2D <- function(funDataObject, bs = "ps", m = NA, k = -1)
{
  if(dimSupp(funDataObject) != 2)
    stop("splines2D is implemented for 2D functional data (images) only.")
  
  N <- nObs(funDataObject)
  
  coord <- expand.grid(x = funDataObject@argvals[[1]], y = funDataObject@argvals[[2]])
  
  # spline design matrix via gam
  g <- mgcv::gam(as.vector(funDataObject@X[1,,]) ~ te(coord$x, coord$y, bs = bs, m = m, k = k), data = coord, fit = FALSE)
  desMat <- g$X
  k <- sapply(g$smooth[[1]]$margin, function(l){l$bs.dim})
  m <- lapply(g$smooth[[1]]$margin, function(l){l$p.order})
  
  # weights via lm -> no penalization
  scores <- t(apply(funDataObject@X, 1, function(f, dM){stats::lm(as.vector(f) ~ dM - 1)$coef}, dM = desMat))
  
  # extract basis functions (in the correct dimensions)
  B <- aperm(array(desMat, c(funData::nObsPoints(funDataObject), ncol(scores))), c(3,1,2))
  
  return(list(scores = scores,
              B = calcBasisIntegrals(B, 2, funDataObject@argvals),
              ortho = FALSE,
              functions = NULL,
              settings = list(bs = bs, k = k, m = m)
  ))
}

#' @rdname splineBasis2D
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#' @importFrom mgcv bam model.matrix.gam te
#' 
#' @export splineBasis2Dpen
splineBasis2Dpen <- function(funDataObject, bs = "ps", m = NA, k = -1, parallel = FALSE)
{
  if(dimSupp(funDataObject) != 2)
    stop("splines2Dpen is implemented for 2D functional data (images) only.")
  
  N <- nObs(funDataObject)
  
  coord <- expand.grid(x = funDataObject@argvals[[1]], y = funDataObject@argvals[[2]])
  
  if(parallel)
  {
    scores <- foreach::foreach(i = 1:(N-1), .combine = "rbind")%dopar%{
      g <- mgcv::bam(as.vector(funDataObject@X[i, , ]) ~ te(coord$x, coord$y, bs = bs, m = m, k = k), data = coord, method = "REML")
      g$coef
    }
  }
  else
  {
    scores <- foreach::foreach(i = 1:(N-1), .combine = "rbind")%do%{
      g <- mgcv::bam(as.vector(funDataObject@X[i, , ]) ~ te(coord$x, coord$y, bs = bs, m = m, k = k), data = coord, method = "REML")
      g$coef
    }
  }
  
  # fit the last one extra in order to extract model matrix
  g <- mgcv::bam(as.vector(funDataObject@X[N, , ]) ~ te(coord$x, coord$y, bs = bs, m = m, k = k), data = coord, method = "REML")
  k <- sapply(g$smooth[[1]]$margin, function(l){l$bs.dim})
  m <- lapply(g$smooth[[1]]$margin, function(l){l$p.order})
  
  scores <- rbind(scores, g$coef)
  
  # extract basis functions (in the correct dimensions)
  B <- aperm(array(mgcv::model.matrix.gam(g), c(nObsPoints(funDataObject), ncol(scores))), c(3,1,2))
  
  return(list(scores = scores,
              B = calcBasisIntegrals(B, 2, funDataObject@argvals),
              ortho = FALSE,
              functions = NULL,
              settings = list(bs = bs, k = k, m = m)
  ))
}


#' Calculate a cosine basis representation for functional data on two- or 
#' three-dimensional domains
#' 
#' These functions calculate  a tensor cosine basis representation for 
#' functional data on two- or three-dimensional domains based on a discrete 
#' cosine transformation (DCT) using the C-library \code{fftw3} 
#' (\url{http://www.fftw.org/}). Coefficients under a given thershold are set to
#' 0 to reduce complexity and for denoising.
#' 
#' Given the (discretized) observed functions \eqn{X_i}, the function 
#' \code{dctBasis2D} calculates a basis representation \deqn{X_i(s,t) = \sum_{m 
#' = 0}^{K_1-1} \sum_{n = 0}^{K_2-1} \theta_{mn} f_{mn}(s,t)} of a 
#' two-dimensional function \eqn{X_i(s,t)} in terms of (orthogonal) tensor 
#' cosine basis functions \deqn{f_{mn}(s,t) = c_m c_n \cos(ms) \cos(nt), \quad 
#' (s,t) \in \mathcal{T}}{f_{mn}(s,t) = c_m c_n \cos(ms) \cos(nt), \quad (s,t) 
#' \in \calT} with \eqn{c_m = \frac{1}{\sqrt{\pi}}} for \eqn{m=0} and \eqn{c_m =
#' \sqrt{\frac{2}{\pi}}} for \eqn{m=1,2,\ldots} based on a discrete cosine 
#' transform (DCT).
#' 
#' If not thresholded (\code{qThresh = 0}), the function returns all non-zero 
#' coefficients \eqn{\theta_{mn}} in the basis representation in a 
#' \code{\link[Matrix]{sparseMatrix}} (package \pkg{Matrix}) called
#' \code{scores}. Otherwise, coefficients with \deqn{|\theta_{mn}| <= q } are
#' set to zero, where \eqn{q} is the \code{qThresh}-quantile of
#' \eqn{|\theta_{mn}|}.
#' 
#' For functions \eqn{X_i(s,t,u)} on three-dimensional domains, the function 
#' \code{dctBasis3D} calculates a basis representation \deqn{X_i(s,t,u) = 
#' \sum_{m = 0}^{K_1-1} \sum_{n = 0}^{K_2-1} \sum_{k = 0}^{K_3-1} \theta_{mnk} 
#' f_{mnk}(s,t,u)} in terms of (orthogonal) tensor cosine basis functions 
#' \deqn{f_{mnk}(s,t,u) = c_m c_n c_k \cos(ms) \cos(nt) \cos(ku), \quad (s,t,u) 
#' \in \mathcal{T}}{f_{mnk}(s,t,u) = c_m c_n c_k \cos(ms) \cos(nt) \cos(ku), 
#' \quad (s,t,u) \in \calT} again with \eqn{c_m = \frac{1}{\sqrt{pi}}} for 
#' \eqn{m=0} and \eqn{c_m = \sqrt{\frac{2}{pi}}} for \eqn{m=1,2,\ldots} based on
#' a discrete cosine transform (DCT). The thresholding works analogous as for 
#' the two-dimensional case.
#' 
#' @section Warning: If the C-library \code{fftw3} is not available when the 
#'   package \code{MFPCA} is installed, this function is disabled an will throw 
#'   an error. For full functionality install the C-library \code{fftw3} from 
#'   \url{http://www.fftw.org/} and reinstall \code{MFPCA}.
#'   
#' @param funDataObject An object of class \code{\link[funData]{funData}} 
#'   containing the observed functional data samples and for which the basis 
#'   representation is calculated.
#' @param qThresh A numeric with value in \eqn{[0,1]}, giving the quantile for 
#'   thresholding the coefficients. See Details.
#' @param parallel Logical. If \code{TRUE}, the coefficients for the basis 
#'   functions are calculated in parallel. The implementation is based on the 
#'   \code{\link[foreach]{foreach}} function and requires a parallel backend 
#'   that must be registered before; see \code{\link[foreach]{foreach}} for 
#'   details. Defaults to \code{FALSE}.
#'   
#' @return \item{scores}{A \code{\link[Matrix]{sparseMatrix}} of scores 
#'   (coefficients) with dimension \code{N x K}, reflecting the weights 
#'   \eqn{\theta_{mn}} (\eqn{\theta_{mnk}}) for each basis function in each
#'   observation, where \code{K} is the total number of basis functions used.}
#'   \item{B}{A diagonal matrix, giving the norms of the different basis
#'   functions used (as they are orthogonal).} \item{ortho}{Logical, set to
#'   \code{FALSE}, as basis functions are orthogonal, but in genereal not
#'   orthonormal.} \item{functions}{\code{NULL}, as basis functions are known.}
#'   
#' @seealso \code{\link{univDecomp}}, \code{\link{dct2D}}, \code{\link{dct3D}}
#'   
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#' @importFrom Matrix sparseMatrix
#'   
#' @export dctBasis2D
#'   
#' @examples
#' # Simulate data with 10 observations on two-dimensional domain (images)
#' x1 <- seq(0, 1, length.out = 50)
#' x2 <- seq(-1, 1, length.out = 75)
#' f2 <- funData(argvals = list(x1, x2),
#'               X = aperm(replicate(10, x1 %o% cos(pi*x2) + 
#'                                   matrix(rnorm(50*75, sd = 0.1), nrow = 50)),
#'                        c(3,1,2)))
#' 
#' # Calculate basis functions: This will throw an error if fftw3 is not installed.           
#' \dontrun{
#' dct2D <- dctBasis2D(f2, qThresh = 0.95)
#' 
#' # verify that scores are saved in a sparse matrix
#' dct2D$scores[,1:25] # the first 25 scores for each observation
#' }
dctBasis2D <- function(funDataObject, qThresh, parallel = FALSE)
{
  if(dimSupp(funDataObject) != 2)
    stop("dctBasis2D can handle only functional data on two-dimensional domains.")
  
  if(parallel)
    res <- foreach::foreach(i = 1:nObs(funDataObject), .combine = "rbind") %dopar% {
      dct <- dct2D(funDataObject@X[i,,], qThresh)
      
      data.frame(i = rep(i, length(dct$ind)), j = dct$ind, x = dct$val)
    }
  else
    res <- foreach::foreach(i = 1:nObs(funDataObject), .combine = "rbind") %do% {
      dct <- dct2D(funDataObject@X[i,,], qThresh)
      
      data.frame(i = rep(i, length(dct$ind)), j = dct$ind, x = dct$val)
    }
  
  return(list(scores = sparseMatrix(i = res$i, j = res$j, x = res$x),
              B = Matrix::Diagonal(n = max(res$j), x = prod(sapply(funDataObject@argvals, function(l){diff(range(l))}))/pi^2),
              ortho = FALSE,
              functions = NULL
  ))
}


#' Calculate and threshold DCT for an image
#'
#' This function calculates the (orthonormal) discrete cosine transformation for
#' an image and returns thresholded DCT coefficients using the C-library
#' \code{fftw3} (see \url{http://www.fftw.org/}).
#'
#' @section Warning: If the C-library \code{fftw3} is not available when the
#'   package \code{MFPCA} is installed, this function is disabled an will throw
#'   an error. For full functionality install the C-library \code{fftw3} from
#'   \url{http://www.fftw.org/} and reinstall \code{MFPCA}.
#'
#' @param image An image (a 2D matrix with real values).
#' @param qThresh A numeric with value in \eqn{[0,1]}, giving the quantile for
#'   thresholding the coefficients. See \code{\link{dctBasis2D}} for details.
#'
#' @return \item{ind}{An integer vector, containing the indices of
#'   non-thresholded (hence non-zero) coefficients.} \item{val}{A numeric
#'   vector, giving the values of the corresponding coefficients.}
#'
#' @seealso \code{\link{dctBasis2D}}
#' 
#' @importFrom stats quantile
#'
#' @useDynLib MFPCA calcCoefs
#'
#' @keywords internal
dct2D <- function(image, qThresh)
{
  if(length(dim(image)) != 2)
    stop("dct2D can handle only 2D images")
  
  res <- .C("calcCoefs", M = as.integer(nrow(image)), N = as.integer(ncol(image)),
            image = as.numeric(image), coefs = as.numeric(image*0))$coefs
  
  ind <- which(abs(res) > stats::quantile(abs(res), qThresh))
  
  return(list(ind = ind, val = res[ind]))
}


#' @rdname dctBasis2D
#' 
#' @export dctBasis3D
dctBasis3D <- function(funDataObject, qThresh, parallel = FALSE)
{
  if(dimSupp(funDataObject) != 3)
    stop("dctBasis2D can handle only functional data on three-dimensional domains.")
  
  if(parallel)
    res <- foreach::foreach(i = 1:nObs(funDataObject), .combine = "rbind") %dopar% {
      dct <- dct3D(funDataObject@X[i,,,], qThresh)
      
      data.frame(i = rep(i, length(dct$ind)), j = dct$ind, x = dct$val)
    }
  else
    res <- foreach::foreach(i = 1:nObs(funDataObject), .combine = "rbind") %do% {
      dct <- dct3D(funDataObject@X[i,,,], qThresh)
      
      data.frame(i = rep(i, length(dct$ind)), j = dct$ind, x = dct$val)
    }
  
  return(list(scores = sparseMatrix(i = res$i, j = res$j, x = res$x),
              B = Matrix::Diagonal(n = max(res$j), x = prod(sapply(funDataObject@argvals, function(l){diff(range(l))}))/pi^3),
              ortho = FALSE,
              functions = NULL
  ))
}


#' Calculate and threshold DCT for an 3D image
#'
#' This function calculates the (orthonormal) discrete cosine transformation for
#' a 3D image and returns thresholded DCT coefficients using the C-library
#' \code{fftw3} (see \url{http://www.fftw.org/}).
#'
#' @section Warning: If the C-library \code{fftw3} is not available when the
#'   package \code{MFPCA} is installed, this function is disabled an will throw
#'   an error. For full functionality install the C-library \code{fftw3} from
#'   \url{http://www.fftw.org/} and reinstall \code{MFPCA}.
#'
#' @param image A 3D image (a 3D array with real values).
#' @param qThresh A numeric with value in \eqn{[0,1]}, giving the quantile for
#'   thresholding the coefficients. See \code{\link{dctBasis3D}} for details.
#'
#' @return \item{ind}{An integer vector, containing the indices of
#'   non-thresholded (hence non-zero) coefficients.} \item{val}{A numeric
#'   vector, giving the values of the corresponding coefficients.}
#'
#' @seealso \code{\link{dctBasis3D}}
#' 
#' @importFrom stats quantile
#'
#' @useDynLib MFPCA calcCoefs3D
#'
#' @keywords internal
dct3D <- function(image, qThresh)
{
  if(length(dim(image)) != 3)
    stop("dct3D can handle only 3D images")
  
  res <- .C("calcCoefs3D", dim = as.integer(dim(image)),
            image = as.numeric(image), coefs = as.numeric(image*0))$coefs
  
  ind <- which(abs(res) > stats::quantile(abs(res), qThresh))
  
  return(list(ind = ind, val = res[ind]))
}


# clean up after unloading
.onUnload <- function (libpath) {
  library.dynam.unload("MFPCA", libpath)
}
