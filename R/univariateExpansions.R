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
#' @seealso \link{MFPCA}, \link{fpcaFunction}, \link{splineFunction1D},
#'   \link{splineFunction2D}, \link{splineFunction2Dpen}, \link{dctFunction2D},
#'   \link{defaultFunction}
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
#' This function calculates a linear combination of arbitrary basis funtions on
#' domains with arbitrary dimension.
#'
#' @param scores A matrix of dimension \eqn{N x K}, representing the \eqn{K}
#'   scores (coefficients) for each observation \eqn{i = 1, \ldots, N}.
#' @param argvals A list containing a vector of x-values, representing a  \eqn{M_1
#'   x ...x M_p} dimensional interval.
#' @param functions A \code{funData} object, representing \eqn{K} basis
#'   functions on a \eqn{M_1 x ...x M_p} dimensional domain.
#'
#' @return An object of class \code{funData} with \eqn{N} observations on
#'   \code{argvals}, corresponding to the linear combination of the basis
#'   functions.
#'
#' @seealso univExpansion
defaultFunction <- function(scores, argvals, functions)
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
#' @param scores A matrix of dimension \eqn{N x K}, representing the \eqn{K}
#'   scores (coefficients) for each observation \eqn{i = 1, \ldots, N}.
#' @param argvals A list containing a vector of x-values.
#' @param functions A \code{funData} object, representing the FPCA basis.
#'
#' @return An object of class \code{funData} with \eqn{N} observations on
#'   \code{argvals}, corresponding to the linear combination of the functional
#'   principal components.
#'
#' @seealso univExpansion
fpcaFunction <- function(scores, argvals, functions)
{
  return(funData(argvals, scores %*% functions@X))
}

#' Calculate a linear combination of a UMPCA basis on two-dimensional domains
#' 
#' This function calculates a linear combination of basis functions based on
#' uncorrelated multilinear principal component analysis (UMPCA) for data on
#' two-dimensional domains.
#' 
#' @param scores A matrix of dimension \eqn{N x K}, representing the \eqn{K} 
#'   scores (coefficients) for each observation \eqn{i = 1, \ldots, N}.
#' @param argvals A list containing a vector of x-values.
#' @param functions A \code{funData} object, representing the UMPCA basis.
#'   
#' @return An object of class \code{funData} with \eqn{N} observations on 
#'   \code{argvals}, corresponding to the linear combination of the functional 
#'   principal components.
#'   
#' @seealso univExpansion
umpcaFunction <- function(scores, argvals, functions)
{
  if(dimSupp(functions) != 2)
    stop("UMPCA option for univExpansion is implemented for 2D data (images) only!")
  
  N <- nrow(scores) # number of observations
  recons <- array(0, c(N, nObsPoints(functions)))
  
  for(i in 1:N)
    recons[i,,] <- ttv(functions@X, list(scores[i,]), dim = 1)
  
  reconsFunctions <- funData(argvals = functions@argvals, X = recons)
}

#' Calculate linear combinations of spline basis functions on one-dimensional
#' domains
#'
#' Given scores (coefficients), this function calculates a linear combination of
#' one-dimensional spline basis functions on one-dimensional domains based on
#' the \link[mgcv]{gam} function in the \pkg{mgcv} package.
#'
#' @param scores A matrix of dimension \eqn{N x K}, representing the \eqn{K}
#'   scores (coefficients) for each observation \eqn{i = 1, \ldots, N}.
#' @param argvals A list containing a vector of x-values.
#' @param bs A character string, specifying the type of basis functions to be
#'   used. Please refer to
#'   \code{\link[mgcv]{smooth.terms}} for a list of possible basis functions.
#' @param m A numeric, the order of the spline basis. See  \code{\link[mgcv]{s}} for details.
#' @param k A numeric, the number of basis functions used. See
#'   \code{\link[mgcv]{s}} for details.
#'
#' @return An object of class \code{funData} with \eqn{N} observations on
#'   \code{argvals}, corresponding to the linear combination of spline basis
#'   functions.
#'
#' @seealso univExpansion
#'
#' @importFrom mgcv gam
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
#' Given scores (coefficients), this function calculates a linear combination of
#' two-dimensional spline tensor basis functions on two-dimensional domains
#' based on the \link[mgcv]{gam} function in the \pkg{mgcv} package. If the
#' scores have been calculated based on a penalized tensor spline basis, use
#' \code{splineFunction2Dpen} instead (which runs \link[mgcv]{bam} instead of
#' \link[mgcv]{gam}).
#'
#' @param scores A matrix of dimension \eqn{N x K}, representing the \eqn{K}
#'   scores (coefficients) for each observation \eqn{i = 1, \ldots, N}.
#' @param argvals A list containing a two numeric vectors, corresponding to the x-
#'   and y-values.
#' @param bs An vector of character strings (or a single character), the type of
#'   basis functions to be used. Please refer to \code{\link[mgcv]{te}} for a
#'   list of possible basis functions.
#' @param m A numeric vector (or a single number), the order of the spline
#'   basis. See \code{\link[mgcv]{s}} for details.
#' @param k A numeric vector (or a single number), the number of basis functions
#'   used.  See  \code{\link[mgcv]{s}} for details.
#'
#' @return An object of class \code{funData} with \eqn{N} observations on the
#'   two-dimensional domain specified by \code{argvals}, corresponding to the
#'   linear combination of spline basis functions.
#'
#' @seealso univExpansion
#'
#' @importFrom mgcv gam
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


#' Calculate linear combinations of orthonormal cosine basis functions on
#' two-dimensional domains
#'
#' Given scores (coefficients), this function calculates a linear combination of
#' two-dimensional cosine tensor basis functions on two-dimensional domains
#' using the C-library \code{fftw3} (see \url{http://www.fftw.org/}).
#'
#' @section Warning: If the C-library \code{fftw3} is not available when the
#'   package \code{MFPCA} is installed, this function is disabled an will throw
#'   an error. For full functionality install the C-library \code{fftw3} from
#'   \url{http://www.fftw.org/} and reinstall \code{MFPCA}.
#'
#' @param  scores A sparse matrix of dimension \eqn{N x L}, representing the
#'   \eqn{L} scores (coefficients) for each observation \eqn{i = 1, \ldots, N}.
#' @param argvals A list containing a two numeric vectors, corresponding to the x-
#'   and y-values.
#' @param parallel Logical. If \code{TRUE}, the coefficients for the basis
#'   functions are calculated in parallel. The implementation is based on the
#'   \code{\link[foreach]{foreach}} function and requires a parallel backend
#'   that must be registered before. See \code{\link[foreach]{foreach}} for
#'   details.
#'
#' @return An object of class \code{funData} with \eqn{N} observations on the
#'   two-dimensional domain specified by \code{argvals}, corresponding to the
#'   linear combination of orthonormal cosine basis functions.
#'
#' @seealso univExpansion
#'
#' @importFrom abind abind
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
#' @seealso dctBasis2D
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


#' Calculate linear combinations of orthonormal cosine basis functions on
#' three-dimensional domains
#'
#' Given scores (coefficients), this function calculates a linear combination of
#' threse-dimensional cosine tensor basis functions on three-dimensional domains
#' using the C-library \code{fftw3} (see \url{http://www.fftw.org/}).
#'
#' @section Warning: If the C-library \code{fftw3} is not available when the
#'   package \code{MFPCA} is installed, this function is disabled an will throw
#'   an error. For full functionality install the C-library \code{fftw3} from
#'   \url{http://www.fftw.org/} and reinstall \code{MFPCA}.
#'
#' @param  scores A sparse matrix of dimension \eqn{N x L}, representing the
#'   \eqn{L} scores (coefficients) for each observation \eqn{i = 1, \ldots, N}.
#' @param argvals A list containing a two numeric vectors, corresponding to the x-
#'   and y-values.
#' @param parallel Logical. If \code{TRUE}, the coefficients for the basis
#'   functions are calculated in parallel. The implementation is based on the
#'   \code{\link[foreach]{foreach}} function and requires a parallel backend
#'   that must be registered before. See \code{\link[foreach]{foreach}} for
#'   details.
#'
#' @return An object of class \code{funData} with \eqn{N} observations on the
#'   two-dimensional domain specified by \code{argvals}, corresponding to the
#'   linear combination of orthonormal cosine basis functions.
#'
#' @seealso univExpansion
#'
#' @importFrom abind abind
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
#' @seealso dctBasis3D
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
