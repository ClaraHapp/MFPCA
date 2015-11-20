#' Univariate Basis Decomposition
#'
#' This function calculates a univariate basis decomposition for a (univariate)
#' functional data object \code{data}.
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
#' @param data A \code{funData} object, representing the (univariate) functional
#'   data samples.
#' @param params A list containing the parameters for the particular basis to
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
#' @seealso \link{MFPCA}, \link{fpcaBasis}, \link{splineBasis1D},
#'   \link{splineBasis1Dpen}, \link{splineBasis2D}, \link{splineBasis2Dpen},
#'   \link{dctBasis2D}
univDecomp <- function(type, data, params)
{
  if(is.null(params))
    params <- list() # create empty list

  params$funDataObject <- data

  res <- switch(type,
                "uFPCA" = do.call(fpcaBasis, params),
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
#' This function calculates functional principal component basis representation
#' for functional data on one-dimensional domains. The FPCA is calculated via
#' the \link{PACE} function, which is built on \link[refund]{fpca.sc} in the
#' \pkg{refund} package.
#'
#' @param funDataObject An object of class \code{\link[funData]{funData}}
#'   containing the observed functional data samples and for which the FPCA is
#'   to be calculated.
#' @param nbasis An integer, representing the number of  B-spline basis
#'   functions used for estimation of the mean function and bivariate smoothing
#'   of the covariance surface. Defaults to 10 (cf.
#'   \code{\link[refund]{fpca.sc}}).
#' @param pve A numeric value between 0 and 1, the proportion of variance
#'   explained: used to choose the number of principal components. Defaults to
#'   0.99 (cf. \code{\link[refund]{fpca.sc}}).
#' @param npc An integer, giving a prespecified value for the number of
#'   principal components. Defaults to \code{NULL}. If given, this overrides
#'   \code{pve} (cf. \code{\link[refund]{fpca.sc}}).
#' @param makePD Logical: should positive definiteness be enforced for the
#'   covariance surface estimate? Defaults to \code{FALSE} (cf.
#'   \code{\link[refund]{fpca.sc}}).
#'
#' @return \item{scores}{A matrix of scores (coefficients) with dimension
#'   \code{N x k}, reflecting the weights for principal component in each
#'   observation.} \item{ortho}{Logical, set to \code{TRUE}, as basis functions are
#'   orthonormal.} \item{functions}{A functional data object, representing the
#'   functional principal component basis functions.}
#'
#' @seealso univDecomp
fpcaBasis <- function(funDataObject, nbasis = 10, pve = 0.99, npc = NULL, makePD = FALSE)
{
  FPCA <- PACE(funDataObject, predData = NULL, nbasis, pve, npc, makePD)

  return(list(scores = FPCA$scores,
              ortho = TRUE,
              functions = FPCA$functions
  ))
}

#' Calculate an unpenalized spline basis decomposition for functional data on
#' one-dimensional domains
#'
#' This function calculates an unpenalized spline basis decomposition for
#' functional data on one-dimensional domains based on the \link[mgcv]{gam}
#' function in the \pkg{mgcv} package.
#'
#' @param funDataObject An object of class \code{\link[funData]{funData}}
#'   containing the observed functional data samples and for which the basis
#'   decomposition is calculated.
#' @param bs A character string, specifying the type of basis functions to be
#'   used. Defaults to "ps" (B-spline functions). Please refer to
#'   \code{\link[mgcv]{smooth.terms}} for a list of possible basis functions.
#' @param m A numeric, the order of the spline basis. Defaults to NA, i.e. the
#'   order is chosen automatically. See  \code{\link[mgcv]{s}} for details.
#' @param k A numeric, the number of basis functions used. Defaults to -1, i.e.
#'   the number of basis functions is chosen automatically. See
#'   \code{\link[mgcv]{s}} for details.
#'
#' @return \item{scores}{A matrix of scores (coefficients) with dimension
#'   \code{N x k}, reflecting the weights for each basis function in each
#'   observation.} \item{B}{A matrix containing the scalar product of all pairs
#'   of basis functions.} \item{ortho}{Logical, set to \code{FALSE}, as basis
#'   functions are not orthonormal.} \item{functions}{\code{NULL}, as basis
#'   functions are known} \item{settings}{A list with entries \code{bs},
#'   \code{m} and \code{k}, giving the actual parameters used for generating the
#'   spline basis functions.}
#'
#' @seealso univDecomp
#'
#' @importFrom mgcv gam
splineBasis1D <- function(funDataObject, bs = "ps", m = NA, k = -1)
{
  N <- nObs(funDataObject)

  x <- funDataObject@xVal[[1]]

  # spline design matrix via gam
  g <- mgcv::gam(funDataObject@X[1, ] ~ s(x, bs = bs, m = m, k = k), fit = FALSE)
  desMat <- g$X
  k <- g$smooth[[1]]$bs.dim
  m <- g$smooth[[1]]$p.order

  # weights via lm -> no penalization
  scores <- t(apply(funDataObject@X, 1, function(f, dM){lm(f ~ dM - 1)$coef}, dM = desMat)) # design matrix already includes intercept!

  return(list(scores = scores,
              B = .calcBasisIntegrals(t(desMat), 1, funDataObject@xVal),
              ortho = FALSE,
              functions = NULL,
              settings = list(bs = bs, k = k, m = m)
  ))
}

#' Calculate a penalized spline basis representation for functional data on
#' one-dimensional domains
#'
#' This function calculates a penalized spline basis representation for
#' functional data on one-dimensional domains based on the \link[mgcv]{gam}
#' function in the \pkg{mgcv} package.
#'
#' @param funDataObject An object of class \code{\link[funData]{funData}}
#'   containing the observed functional data samples and for which the basis
#'   representation is calculated.
#' @param bs A character string, specifying the type of basis functions to be
#'   used. Defaults to "ps" (P-spline functions). Please refer to
#'   \code{\link[mgcv]{smooth.terms}} for a list of possible basis functions.
#' @param m A numeric, the order of the spline basis. Defaults to NA, i.e. the
#'   order is chosen automatically. See  \code{\link[mgcv]{s}} for details.
#' @param k A numeric, the number of basis functions used. Defaults to -1, i.e.
#'   the number of basis functions is chosen automatically. See
#'   \code{\link[mgcv]{s}} for details.
#' @param parallel Logical. If \code{TRUE}, the coefficients for the basis
#'   functions are calculated in parallel. The implementation is based on the
#'   \code{\link[foreach]{foreach}} function and requires a parallel backend
#'   that must be registered before. See \code{\link[foreach]{foreach}} for
#'   details.
#'
#' @return \item{scores}{A matrix of scores (coefficients) with dimension
#'   \code{N x k}, reflecting the weights for each basis function in each
#'   observation.} \item{B}{A matrix containing the scalar product of all pairs
#'   of basis functions.} \item{ortho}{Logical, set to \code{FALSE}, as basis
#'   functions are not orthonormal.} \item{functions}{\code{NULL}, as basis
#'   functions are known} \item{settings}{A list with entries \code{bs},
#'   \code{m} and \code{k}, giving the actual parameters used for generating the
#'   spline basis functions.}
#'
#' @seealso univDecomp
#'
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#' @importFrom mgcv gam
splineBasis1Dpen <- function(funDataObject, bs = "ps", m = NA, k = -1, parallel = FALSE)
{
  N <- nObs(funDataObject)

  x <- funDataObject@xVal[[1]]

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
              B = .calcBasisIntegrals(t(model.matrix(g)), 1, funDataObject@xVal),
              ortho = FALSE,
              functions = NULL,
              settings = list(bs = bs, k = k, m = m)
  ))
}


#' Calculate an unpenalized spline basis representation for functional data on
#' two-dimensional domains
#'
#' This function calculates an unpenalized tensor product spline basis
#' representation for functional data on two-dimensional domains based on the
#' \link[mgcv]{gam} function in the \pkg{mgcv} package.
#'
#' @param funDataObject An object of class \code{\link[funData]{funData}}
#'   containing the observed functional data samples and for which the basis
#'   representation is calculated.
#' @param bs A vector of character strings (or a single character string),
#'   specifying the type of basis functions to be used. Defaults to "ps"
#'   (P-spline functions). Please refer to \code{\link[mgcv]{te}} for a list of
#'   possible basis functions.
#' @param m A numeric vector (or a single number), the order of the spline
#'   basis. Defaults to NA, i.e. the order is chosen automatically.  See
#'   \code{\link[mgcv]{s}} for details.
#' @param k An numeric vector (or a single number), the number of basis
#'   functions used.  Defaults to -1, i.e. the number of basis functions is
#'   chosen automatically.   See  \code{\link[mgcv]{s}} for details.
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
#' @seealso univDecomp
#'
#' @importFrom mgcv gam
splineBasis2D <- function(funDataObject, bs = "ps", m = NA, k = -1)
{
  N <- nObs(funDataObject)

  coord <- expand.grid(x = funDataObject@xVal[[1]], y = funDataObject@xVal[[2]])

  # spline design matrix via gam
  g <- mgcv::gam(as.vector(funDataObject@X[1,,]) ~ te(coord$x, coord$y, bs = bs, m = m, k = k), data = coord, fit = FALSE)
  desMat <- g$X
  k <- sapply(g$smooth[[1]]$margin, function(l){l$bs.dim})
  m <- lapply(g$smooth[[1]]$margin, function(l){l$p.order})

  # weights via lm -> no penalization
  scores <- t(apply(funDataObject@X, 1, function(f, dM){lm(as.vector(f) ~ dM - 1)$coef}, dM = desMat))

  # extract basis functions (in the correct dimensions)
  B <- aperm(array(desMat, c(funData::nObsPoints(funDataObject), ncol(scores))), c(3,1,2))

  return(list(scores = scores,
              B = .calcBasisIntegrals(B, 2, funDataObject@xVal),
              ortho = FALSE,
              functions = NULL,
              settings = list(bs = bs, k = k, m = m)
  ))
}

#' Calculate a penalized spline basis representation for functional data on
#' two-dimensional domains
#'
#' This function calculates a penalized tensor product spline basis
#' representation for functional data on two-dimensional domains based on the
#' \link[mgcv]{bam} function in the \pkg{mgcv} package (for large GAMs).
#'
#' @param funDataObject An object of class \code{\link[funData]{funData}}
#'   containing the observed functional data samples and for which the basis
#'   representation is calculated.
#' @param bs A vector of character strings (or a single character string),
#'   specifying the type of basis functions to be used. Defaults to "ps"
#'   (P-spline functions). Please refer to \code{\link[mgcv]{te}} for a list of
#'   possible basis functions.
#' @param m A numeric vector (or a single number), the order of the spline
#'   basis. Defaults to NA, i.e. the order is chosen automatically.  See
#'   \code{\link[mgcv]{s}} for details.
#' @param k An numeric vector (or a single number), the number of basis
#'   functions used.  Defaults to -1, i.e. the number of basis functions is
#'   chosen automatically.   See  \code{\link[mgcv]{s}} for details.
#' @param parallel Logical. If \code{TRUE}, the coefficients for the basis
#'   functions are calculated in parallel. The implementation is based on the
#'   \code{\link[foreach]{foreach}} function and requires a parallel backend
#'   that must be registered before. See \code{\link[foreach]{foreach}} for
#'   details.
#'
#' @return \item{scores}{A matrix of scores (coefficients) with dimension
#'   \code{N x K}, reflecting the weights for each basis function in each
#'   observation, where \code{K} is the total number of basis functions used.}
#'   \item{B}{A matrix containing the scalar product of all pairs of basis
#'   functions.} \item{ortho}{Logical, set to \code{FALSE}, as basis functions
#'   are not orthonormal.} \item{functions}{\code{NULL}, as basis functions are
#'   known.} \item{settings}{A list with entries \code{bs},
#'   \code{m} and \code{k}, giving the actual parameters used for generating the
#'   spline basis functions.}
#'
#' @seealso univDecomp
#'
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#' @importFrom mgcv bam
splineBasis2Dpen <- function(funDataObject, bs = "ps", m = NA, k = -1, parallel = FALSE)
{
  N <- nObs(funDataObject)

  coord <- expand.grid(x = funDataObject@xVal[[1]], y = funDataObject@xVal[[2]])

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
  B <- aperm(array(model.matrix(g), c(nObsPoints(funDataObject), ncol(scores))), c(3,1,2))

  return(list(scores = scores,
              B = .calcBasisIntegrals(B, 2, funDataObject@xVal),
              ortho = FALSE,
              functions = NULL,
              settings = list(bs = bs, k = k, m = m)
  ))
}


#' Calculate a cosine basis representation for functional data on
#' two-dimensional domains
#'
#' This function calculates a tensor cosine basis representation for functional
#' data on two-dimensional domains based on a discrete cosine transformation
#' (DCT) using the C-library \code{fftw3} \url{http://www.fftw.org/}.
#'
#' Given the (discretized) observed functions \eqn{X_i}, this function
#' calculates a basis representation \deqn{X_i(s,t) = \sum_{m = 0}^M \sum_{n =
#' 0}^N \theta_{mn} f_{mn}(s,t)} in terms of (orthogonal) tensor cosine basis
#' functions \deqn{f_{mn}(s,t) = c_m c_n \cos(ms) \cos(nt), \quad (s,t) \in
#' \mathcal{T}}{f_{mn}(s,t) = c_m c_n \cos(ms) \cos(nt), \quad (s,t) \in
#' \calT} with \eqn{c_m = \frac{1}{\sqrt{pi}}} for \eqn{m=0} and \eqn{c_m =
#' \sqrt{\frac{2}{pi}}} for \eqn{m=1,2,\ldots} based on a discrete cosine
#' transform (DCT).
#'
#' If not thresholded (\code{qThresh = 0}), the function returns all non-zero
#' coefficients \eqn{\theta_{mn}} in the basis representation in a sparse matrix
#' \code{scores}. Otherwise, coefficients with \deqn{|\theta_{mn}| <= q } are
#' set to zero, where \eqn{q} is the \code{qThresh}-quantile of
#' \eqn{|\theta_{mn}|}.
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
#'   that must be registered before. See \code{\link[foreach]{foreach}} for
#'   details.
#'
#' @return \item{scores}{A sparse matrix of scores (coefficients) with dimension
#'   \code{N x L}, reflecting the weights \eqn{\theta_{mn}} for each basis
#'   function in each observation, where \code{L} is the total number of basis
#'   functions used.} \item{B}{A diagonal matrix, giving the norms of the
#'   different basis functions used (as they are orthogonal).}
#'   \item{ortho}{Logical, set to \code{FALSE}, as basis functions are
#'   orthogonal, but in genereal not orthonormal.} \item{functions}{\code{NULL},
#'   as basis functions are known.}
#'
#' @seealso univDecomp
#'
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#' @importFrom Matrix sparseMatrix
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
              B = Matrix::Diagonal(n = max(res$j), x = prod(sapply(funDataObject@xVal, function(l){diff(range(l))}))/pi^2),
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
#'   thresholding the coefficients. See \link{dctBasis2D} for details.
#'
#' @return \item{ind}{An integer vector, containing the indices of
#'   non-thresholded (hence non-zero) coefficients.} \item{val}{A numeric
#'   vector, giving the values of the corresponding coefficients.}
#'
#' @seealso dctBasis2D
#'
#' @useDynLib MFPCA calcCoefs
#'
#' @keywords internal
dct2D <- function(image, qThresh)
{
  res <- .C("calcCoefs", M = as.integer(nrow(image)), N = as.integer(ncol(image)),
            image = as.numeric(image), coefs = as.numeric(image*0))$coefs

  ind <- which(abs(res) > quantile(abs(res), qThresh))

  return(list(ind = ind, val = res[ind]))
}


#' Calculate a cosine basis representation for functional data on
#' three-dimensional domains
#'
#' This function calculates a tensor cosine basis representation for functional
#' data on three-dimensional domains based on a discrete cosine transformation
#' (DCT) using the C-library \code{fftw3} \url{http://www.fftw.org/}.
#'
#' Given the (discretized) observed functions \eqn{X_i}, this function
#' calculates a basis representation \deqn{X_i(s,t, u) = \sum_{m = 0}^{M-1}
#' \sum_{n = 0}^{N-1} \sum_{k = 0}^{K-1} \theta_{mnk}  f_{mnk}(s,t,u)} in terms
#' of (orthogonal) tensor cosine basis functions \deqn{f_{mnk}(s,t,u) = c_m c_n
#' c_k \cos(ms) \cos(nt) \cos(ku), \quad (s,t,u) \in \mathcal{T}}{f_{mnk}(s,t,u) = c_m c_n
#' c_k \cos(ms) \cos(nt) \cos(ku), \quad (s,t,u) \in \calT} with \eqn{c_m =
#' \frac{1}{\sqrt{pi}}} for \eqn{m=0} and \eqn{c_m = \sqrt{\frac{2}{pi}}} for
#' \eqn{m=1,2,\ldots} based on a discrete cosine transform (DCT).
#'
#' If not thresholded (\code{qThresh = 0}), the function returns all non-zero
#' coefficients \eqn{\theta_{mnk}} in the basis representation in a sparse
#' matrix \code{scores}. Otherwise, coefficients with \deqn{|\theta_{mnk}| <= q
#' } are set to zero, where \eqn{q} is the \code{qThresh}-quantile of
#' \eqn{|\theta_{mnk}|}.
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
#'   that must be registered before. See \code{\link[foreach]{foreach}} for
#'   details.
#'
#' @return \item{scores}{A sparse matrix of scores (coefficients) with dimension
#'   \code{N x L}, reflecting the weights \eqn{\theta_{mnk}} for each basis
#'   function in each observation, where \code{L} is the total number of basis
#'   functions used.} \item{B}{A diagonal matrix, giving the norms of the
#'   different basis functions used (as they are orthogonal).}
#'   \item{ortho}{Logical, set to \code{FALSE}, as basis functions are
#'   orthogonal, but in genereal not orthonormal.} \item{functions}{\code{NULL},
#'   as basis functions are known.}
#'
#' @seealso univDecomp
#'
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#' @importFrom Matrix sparseMatrix
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
              B = Matrix::Diagonal(n = max(res$j), x = prod(sapply(funDataObject@xVal, function(l){diff(range(l))}))/pi^3),
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
#'   thresholding the coefficients. See \link{dctBasis3D} for details.
#'
#' @return \item{ind}{An integer vector, containing the indices of
#'   non-thresholded (hence non-zero) coefficients.} \item{val}{A numeric
#'   vector, giving the values of the corresponding coefficients.}
#'
#' @seealso dctBasis3D
#'
#' @useDynLib MFPCA calcCoefs3D
#'
#' @keywords internal
dct3D <- function(image, qThresh)
{
  res <- .C("calcCoefs3D", dim = as.integer(dim(image)),
            image = as.numeric(image), coefs = as.numeric(image*0))$coefs

  ind <- which(abs(res) > quantile(abs(res), qThresh))

  return(list(ind = ind, val = res[ind]))
}


# clean up after unloading
.onUnload <- function (libpath) {
  library.dynam.unload("MFPCA", libpath)
}
