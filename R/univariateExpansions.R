# define global variable i, used by the foreach package and confusing R CMD CHECK
globalVariables('i')

#' @import funData
NULL

#' Calculate a univariate basis expansion
#'
#' This function calculates a univariate basis expansion based on given scores
#' (coefficients) and basis functions.
#'
#' This function calculates functional data \eqn{X_i(t), i= 1 \ldots N} that is
#' represented as a linear combination of basis functions \eqn{b_k(t)}
#' \deqn{X_i(t) = \sum_{k = 1}^K \theta_{ik} b_k(t), i = 1, \ldots, N.} The
#' basis functions may be prespecified (such as spline basis functions or
#' Fourier bases) or can be estimated from observed data (e.g. by functional
#' principal component analysis). If \code{type = "default"} (i.e. a linear
#' combination of arbitrary basis functions is to be calculated), the scores and
#' basis functions must be supplied.
#'
#' @param type A character string, specifying the basis for which the
#'   decomposition is to be calculated.
#' @param argvals A list, representing the domain of the basis functions. See
#'   \linkS4class{funData} for details.
#' @param scores A matrix of scores (coefficients) for each observation based on
#'   the given basis functions.
#' @param functions A functional data object, representing the basis functions.
#'   Can be \code{NULL} if the basis functions are not estimated from observed
#'   data, but have a predefined form. See Details.
#' @param params A list containing the parameters for the particular basis to
#'   use.
#'
#' @return A functional data object, representing the linear combinations of the
#'   basis functions based on the given scores.
#'
#' @seealso \code{\link{MFPCA}}, \code{\link{fpcaFunction}}, \code{\link{splineFunction1D}},
#'   \code{\link{splineFunction2D}}, \code{\link{splineFunction2Dpen}}, \code{\link{dctFunction2D}},
#'  \code{\link{defaultFunction}}
#'
#' @export univExpansion
univExpansion <- function(type, scores, argvals, functions, params = NULL)
{
  params$scores <- scores
  params$functions <- functions

  if(is.numeric(argvals))
    argvals <- list(argvals)

  params$argvals <- argvals

  res <- switch(type,
                "uFPCA" = do.call(fpcaFunction, params),
                "UMPCA" = do.call(umpcaFunction, params),
                "FCP_TPA" = do.call(fcptpaFunction, params),
                "splines1D" = do.call(splineFunction1D, params),
                "splines1Dpen" = do.call(splineFunction1D, params),
                "splines2D" = do.call(splineFunction2D, params),
                "splines2Dpen" = do.call(splineFunction2Dpen, params),
                "DCT2D" = do.call(dctFunction2D, params),
                "DCT3D" = do.call(dctFunction3D, params),
                "default" = do.call(defaultFunction, params),
                stop("Univariate Expansion for 'type' = ", type, " not defined!")
  )

  return(res)
}


#' Calculate a linear combination of arbitrary basis function
#' 
#' This function calculates a linear combination of arbitrary basis functions on 
#' domains with arbitrary dimension.
#' 
#' @param scores A matrix of dimension \code{N x K}, representing the \code{K} 
#'   scores (coefficients) for each of the \code{N} observations.
#' @param argvals A list representing the domain, see \code{\link[funData]{funData}}
#'   for details. Defaults to \code{functions@@argvals}.
#' @param functions A \code{funData} object, representing \code{K} basis 
#'   functions on a domain with arbitrary dimension.
#'   
#' @return An object of class \code{funData} with \code{N} observations on 
#'   \code{argvals}, corresponding to the linear combination of the basis 
#'   functions.
#'   
#' @seealso \code{\link{univExpansion}}
#' 
#' @export defaultFunction
#' 
#' @examples
#' N <- 12 # 12 observations
#' K <- 10 # 10 basis functions
#' 
#' # Use the first 10 Fourier basis functions on [0,1]
#' x <- seq(0,1,0.01)
#' b <- eFun(argvals = x, M = K, type = "Fourier")
#' 
#' # generate N x K score matrix
#' scores <- t(replicate(N, rnorm(K, sd = 1 / (1:10))))
#' 
#' # calculate basis expansion
#' f <- defaultFunction(scores = scores, functions = b) # default value for argvals
#' 
#' oldpar <- par(no.readonly = TRUE)
#' 
#' par(mfrow = c(1,2))
#' plot(b, main = "Basis functions")
#' plot(f, main = "Linear combination")
#' 
#' par(oldpar)
defaultFunction <- function(scores, argvals = functions@argvals, functions)
{
  d <- dim(functions@X)

  # collapse higher-dimensional functions and resize the result
  return(funData(argvals, array(scores %*% array(functions@X, dim = c(d[1], prod(d[-1]))),
                             dim = c(dim(scores)[1],d[-1]))) )
}

#' Calculate a linear combination of a functional principal component basis on 
#' one-dimensional domains
#' 
#' This function calculates  a linear combination of functional principal 
#' component basis functions on one-dimensional domains.
#' 
#' @param scores A matrix of dimension \code{N x K}, representing the \code{K} 
#'   scores (coefficients) for each of the \code{N} observations.
#' @param argvals A list containing a vector of x-values. Defaults to
#'   \code{functions@@argvals}.
#' @param functions A \code{funData} object, representing the FPCA basis.
#'   
#' @return An object of class \code{funData} with \code{N} observations on 
#'   \code{argvals}, corresponding to the linear combination of the functional 
#'   principal components.
#'   
#' @seealso \code{\link{univExpansion}}, \code{\link{fpcaBasis}}
#' 
#' @export fpcaFunction
fpcaFunction <- function(scores, argvals = functions@argvals, functions)
{
  return(funData(argvals, scores %*% functions@X))
}

#' Calculate a linear combination of a UMPCA basis on two-dimensional domains
#' 
#' This function calculates a linear combination of basis functions based on
#' uncorrelated multilinear principal component analysis (UMPCA) for data on
#' two-dimensional domains.
#' 
#' @param scores A matrix of dimension \code{N x K}, representing the \code{K} 
#'   scores (coefficients) for each of the \code{N} observations.
#' @param argvals  A list containing a vector of x-values. Defaults to
#'   \code{functions@@argvals}.
#' @param functions A \code{funData} object, representing the UMPCA basis.
#'   
#' @return An object of class \code{funData} with \code{N} observations on 
#'   \code{argvals}, corresponding to the linear combination of the functional 
#'   principal components.
#'   
#' @seealso \code{\link{univExpansion}}
umpcaFunction <- function(scores, argvals = functions@argvals, functions)
{
  if(dimSupp(functions) != 2)
    stop("UMPCA option for univExpansion is implemented for 2D data (images) only!")
  
  N <- nrow(scores) # number of observations
  recons <- array(0, c(N, nObsPoints(functions)))
  
  for(i in 1:N)
    recons[i,,] <- ttv(functions@X, list(scores[i,]), dim = 1)
  
  reconsFunctions <- funData(argvals = functions@argvals, X = recons)
}


fcptpaFunction <- function(scores, argvals, functions)
{
  return(defaultFunction(scores, argvals, functions))
}

#' Calculate linear combinations of spline basis functions on one-dimensional 
#' domains
#' 
#' Given scores (coefficients), this function calculates a linear combination of
#' spline basis functions on one-dimensional domains based on the 
#' \link[mgcv]{gam} function in the \pkg{mgcv} package.
#' 
#' @param scores A matrix of dimension \code{N x K}, representing the \code{K} 
#'   scores (coefficients) for each of the \code{N} observations.
#' @param argvals A list containing a vector of x-values, on which the functions
#'   should be defined.
#' @param bs A character string, specifying the type of basis functions to be 
#'   used. Please refer to \code{\link[mgcv]{smooth.terms}} for a list of 
#'   possible basis functions.
#' @param m A numeric, the order of the spline basis. See  \code{\link[mgcv]{s}}
#'   for details.
#' @param k A numeric, the number of basis functions used. See 
#'   \code{\link[mgcv]{s}} for details.
#'   
#' @return An object of class \code{funData} with \code{N} observations on 
#'   \code{argvals}, corresponding to the linear combination of spline basis 
#'   functions.
#'   
#' @seealso \code{\link{univExpansion}}, \code{\link{gam}},
#'   \code{\link{splineBasis1D}}
#'   
#' @importFrom mgcv gam
#' 
#' @export splineFunction1D
splineFunction1D <- function(scores, argvals, bs, m, k)
{
  N <- nrow(scores)

  x <- argvals[[1]]

  # spline design matrix via gam
  desMat <- mgcv::gam(rep(0, length(x)) ~ s(x, bs = bs, m = m, k = k), fit = FALSE)$X

  # calculate functions as linear combination of splines
  res <- funData(argvals,
                 scores %*% t(desMat))

  return(res)
}


#' Calculate linear combinations of spline basis functions on two-dimensional 
#' domains
#' 
#' Given scores (coefficients), these functions calculate a linear combination 
#' of spline tensor basis functions on two-dimensional domains based on the 
#' \code{\link[mgcv]{gam}}/\code{\link[mgcv]{bam}} functions in the \pkg{mgcv} package. See 
#' Details.
#' 
#' If the scores have been calculated based on an unpenalized tensor spline 
#' basis, the linear combination is computed based on the
#' \code{\link[mgcv]{gam}} functions ((\code{splineFunction2D})). If the scores
#' were obtained using penalization, the expansion is calculated via 
#' \link[mgcv]{bam} (\code{splineFunction2Dpen}).
#' 
#' @param scores A matrix of dimension \code{N x K}, representing the \code{K} 
#'   scores (coefficients) for each of the \code{N} observations.
#' @param argvals A list containing a two numeric vectors, corresponding to the 
#'   x- and y-values, on which the functions should be defined.
#' @param bs A vector of character strings (or a single character), the type of 
#'   basis functions to be used. Please refer to \code{\link[mgcv]{te}} for a 
#'   list of possible basis functions.
#' @param m A numeric vector (or a single number), the order of the spline 
#'   basis. See \code{\link[mgcv]{s}} for details.
#' @param k A numeric vector (or a single number), the number of basis functions
#'   used.  See  \code{\link[mgcv]{s}} for details.
#'   
#' @return An object of class \code{funData} with \code{N} observations on the 
#'   two-dimensional domain specified by \code{argvals}, corresponding to the 
#'   linear combination of spline basis functions.
#'   
#' @seealso \code{\link{univExpansion}},  \code{\link{gam}}, 
#'   \code{\link{splineBasis2D}}
#'   
#' @importFrom mgcv gam
#'   
#' @export splineFunction2D
splineFunction2D <- function(scores, argvals, bs, m, k)
{
  N <- nrow(scores)

  coord <- expand.grid(x = argvals[[1]], y = argvals[[2]])

  # spline design matrix via gam
  desMat <- mgcv::gam(rep(0, dim(coord)[1]) ~ te(coord$x, coord$y, bs = bs, m = m, k = k), data = coord, fit = FALSE)$X

  # calculate functions as linear combination of splines
  res <- funData(argvals,
                 array(scores %*% t(desMat),
                       dim = c(N, length(argvals[[1]]), length(argvals[[2]]))))

  return(res)
}


#' @rdname splineFunction2D
#'   
#' @importFrom mgcv gam
#'   
#' @export splineFunction2Dpen
splineFunction2Dpen <- function(scores, argvals, bs, m, k)
{
  N <- nrow(scores)

  coord <- expand.grid(x = argvals[[1]], y = argvals[[2]])

  # spline design matrix via gam
  desMat <- mgcv::bam(rep(0, dim(coord)[1]) ~ te(coord$x, coord$y, bs = bs, m = m, k = k), data = coord, fit = FALSE,
                      chunk.size = nrow(coord))$X # use exactly the given grid, no chunks

  # calculate functions as linear combination of splines
  res <- funData(argvals,
                 array(scores %*% t(desMat),
                       dim = c(N, length(argvals[[1]]), length(argvals[[2]]))))

  return(res)
}


#' Calculate linear combinations of orthonormal cosine basis functions on two- 
#' or three-dimensional domains
#' 
#' Given scores (coefficients), these functions calculate a linear combination 
#' of two- or three-dimensional cosine tensor basis functions on two- or 
#' three-dimensional domains using the C-library \code{fftw3} (see 
#' \url{http://www.fftw.org/}).
#' 
#' @section Warning: If the C-library \code{fftw3} is not available when the 
#'   package \code{MFPCA} is installed, the functions are disabled an will throw
#'   an error. For full functionality install the C-library \code{fftw3} from 
#'   \url{http://www.fftw.org/} and reinstall \code{MFPCA}.
#'   
#' @param  scores A sparse matrix of dimension \code{N x L}, representing the 
#'   \code{L} scores (coefficients), where \code{N} is the number of 
#'   observations.
#' @param argvals A list containing two or three numeric vectors, corresponding 
#'   to the domain grid (x and y values for two-dimensional domains; x,y and z 
#'   values fro three-dimensional domains.)
#' @param parallel Logical. If \code{TRUE}, the coefficients for the basis 
#'   functions are calculated in parallel. The implementation is based on the 
#'   \code{\link[foreach]{foreach}} function and requires a parallel backend 
#'   that must be registered before; see \code{\link[foreach]{foreach}} for 
#'   details. Defaults to \code{FALSE}.
#'   
#' @return An object of class \code{funData} with \code{N} observations on the 
#'   two- or threedimensional domain specified by \code{argvals}, corresponding 
#'   to the linear combination of orthonormal cosine basis functions.
#'   
#' @seealso \code{\link{univExpansion}}, \code{\link{idct2D}}, 
#'   \code{\link{idct3D}}, \code{\link{dctBasis2D}}, \code{\link{dctBasis3D}}
#'   
#' @importFrom abind abind
#'   
#' @export dctFunction2D
#'   
#' @examples
#' # generate sparse 10 x 15 score matrix (i.e. 10 observations) with 30 entries
#' scores <- Matrix::sparseMatrix(i = sample(1:10, 30, replace = TRUE), # sample row indices
#'      j = sample(1:15, 30, replace = TRUE, prob = 1/(1:15)), # sample column indices, high indices with low probability
#'      x = rnorm(30)) # sample values
#' scores
#' 
#' # calculate basis expansion on [0,1] x [0,1]
#' f <- dctFunction2D(scores = scores, argvals = list(seq(0,1,0.01), seq(0,1,0.01)))
#' nObs(f) # f has 10 observations
#' plot(f, obs = 1) # plot first observation
#' plot(f, obs = 2) # plot second observation
dctFunction2D <- function(scores, argvals, parallel = FALSE)
{
  # dimension of the image
  dim <-sapply(argvals, length)

  # get indices of sparse matrix
  scores <- as(scores, "dgTMatrix") # uncompressed format

  if(parallel)
    res <- foreach::foreach(i = 0:max(scores@i), .combine = function(x,y){abind(x,y, along = 3)}) %dopar%{
      idct2D(scores@x[scores@i == i], 1 + scores@j[scores@i == i], dim = dim) # 0-indexing!
    }
  else
    res <- foreach::foreach(i = 0:max(scores@i), .combine =function(x,y){abind(x,y, along = 3)}) %do%{
      idct2D(scores@x[scores@i == i], 1 + scores@j[scores@i == i], dim = dim) # 0-indexing!
    }

  return(funData(argvals, X = aperm(res, c(3,1,2))))
}


#' Calculate an inverse DCT for an image
#'
#' This function calculates an inverse (orthonormal) discrete cosine
#' transformation for given coefficients in two dimensions using the C-library
#' \code{fftw3} (see \url{http://www.fftw.org/}). As many coefficients are
#' expected to be zero, the values are given in compressed format (indices and
#' values only of non-zero coefficients).
#'
#' @section Warning: If the C-library \code{fftw3} is not available when the
#'   package \code{MFPCA} is installed, this function is disabled an will throw
#'   an error. For full functionality install the C-library \code{fftw3} from
#'   \url{http://www.fftw.org/} and reinstall \code{MFPCA}.
#'
#' @param scores A numeric vector, containing the non-zero coefficients.
#' @param ind An integer vector, containing the indices of the non-zero
#'   coefficients.
#' @param dim A numeric vector of length 2, giving the resulting image
#'   dimensions.
#'
#' @return A matrix of dimensions \code{dim}, which is a linear combination of
#'   cosine tensor basis functions with the given coefficients.
#'
#' @seealso \code{\link{dctBasis2D}}
#'
#' @useDynLib MFPCA calcImage
#'
#' @keywords internal
idct2D <- function(scores, ind, dim)
{
  if(length(dim) != 2)
    stop("Function idct2D can handle only 2D images.")
  
  if(min(ind) < 1)
    stop("Indices must be positive.")
  
  if(max(ind) > prod(dim))
    stop("Index exceeds image dimensions.")
  
  full <- array(0, dim)
  full[ind] <- scores

  res <- .C("calcImage", M = as.integer(dim[1]), N = as.integer(dim[2]),
            coefs = as.numeric(full), image = as.numeric(full*0))$image

  return(array(res, dim))
}


#' @rdname dctFunction2D
#'   
#' @export dctFunction3D
dctFunction3D <- function(scores, argvals, parallel = FALSE)
{
  # dimension of the image
  dim <-sapply(argvals, length)

  # get indices of sparse matrix
  scores <- as(scores, "dgTMatrix") # uncompressed format

  if(parallel)
    res <- foreach::foreach(i = 0:max(scores@i), .combine = function(x,y){abind(x,y, along = 4)}) %dopar%{
      idct3D(scores@x[scores@i == i], 1 + scores@j[scores@i == i], dim = dim) # 0-indexing!
    }
  else
    res <- foreach::foreach(i = 0:max(scores@i), .combine =function(x,y){abind(x,y, along = 4)}) %do%{
      idct3D(scores@x[scores@i == i], 1 + scores@j[scores@i == i], dim = dim) # 0-indexing!
    }

  return(funData(argvals, X = aperm(res, c(4,1,2,3))))
}


#' Calculate an inverse DCT for a 3D image
#'
#' This function calculates an inverse (orthonormal) discrete cosine
#' transformation for given coefficients in three dimensions using the C-library
#' \code{fftw3} (see \url{http://www.fftw.org/}). As many coefficients are
#' expected to be zero, the values are given in compressed format (indices and
#' values only of non-zero coefficients).
#'
#' @section Warning: If the C-library \code{fftw3} is not available when the
#'   package \code{MFPCA} is installed, this function is disabled an will throw
#'   an error. For full functionality install the C-library \code{fftw3} from
#'   \url{http://www.fftw.org/} and reinstall \code{MFPCA}.
#'
#' @param scores A numeric vector, containing the non-zero coefficients.
#' @param ind An integer vector, containing the indices of the non-zero
#'   coefficients.
#' @param dim A numeric vector of length 3, giving the resulting image
#'   dimensions.
#'
#' @return A matrix of dimensions \code{dim}, which is a linear combination of
#'   cosine tensor basis functions with the given coefficients.
#'
#' @seealso \code{\link{dctBasis3D}}
#'
#' @useDynLib MFPCA calcImage3D
#'
#' @keywords internal
idct3D <- function(scores, ind, dim)
{
  if(length(dim) != 3)
    stop("Function idct3D can handle only 3D images.")
  
  if(min(ind) < 1)
    stop("Indices must be positive.")
  
  if(max(ind) > prod(dim))
    stop("Index exceeds image dimensions.")
  
  full <- array(0, dim)
  full[ind] <- scores

  res <- .C("calcImage3D", dim = as.integer(dim),
            coefs = as.numeric(full), image = as.numeric(full*0))$image

  return(array(res, dim))
}
