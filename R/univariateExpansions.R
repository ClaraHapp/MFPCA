# define global variable i, used by the foreach package and confusing R CMD CHECK
globalVariables('i')

#' @import funData
NULL

#' Calculate a univariate spline basis expansion for a functional data object
#'
#' This function calculates a univariate spline basis expansion for a functional
#' data object defined on a one- or two-dimensional domain. All basis functions
#' implemented for the functions \code{\link[mgcv]{gam}} /
#' \code{\link[mgcv]{bam}} in the \pkg{mgvc} package can be used. For functions
#' with two-dimensional domain, spline basis tensor products are used. If data
#' should be pre-smoothed, e.g. due to considerable measurement error, the
#' expansion is fit using an appropriate penalty term.
#'
#' @param funDataObject An object of class \code{\link[funData]{funData}}
#'   containing the observed functional data, for which the bases representation
#'   is calculated.
#' @param bs A character string, the type of basis functions to be used. Defaults to p-Splines
#'   (\code{"ps"}) or tensor products of p-splines in case of two-dimensional support. Please
#'   refer to \code{\link[mgcv]{smooth.terms}} for a list of possible basis
#'   functions.
#' @param m An array (or a single character), the order of the spline basis.
#'   Defaults to \code{2}, as in \code{\link[mgcv]{gam}} /
#'   \code{\link[mgcv]{bam}}. See \code{\link[mgcv]{s}} for details.
#' @param k An array (or a single character), the number of basis functions
#'   used.  Defaults to \code{-1}, as in \code{\link[mgcv]{gam}} /
#'   \code{\link[mgcv]{bam}}. See \code{\link[mgcv]{s}} for details.
#' @param pen Logical. If \code{TRUE}, the basis expansion is fit using an
#'   appropriate penalization term. Defaults to \code{FALSE}.
#'
#' @return \item{scores}{A matrix of weights with dimension \code{N x k},
#' reflecting the weights for each basis function in each observation.}
#' \item{functions}{A \code{\link[funData]{funData}} object representing the
#' basis functions.} \item{basisLong}{A matrix representing the basis functions
#' in vectorized form (only if \code{funDataobject} is defined on a
#' two-dimensional domain).}
#'
#' @seealso \code{\link[mgcv]{gam}}, \code{\link[mgcv]{s}},
#'   \code{\link[mgcv]{te}}, \code{\link{MFPCA}}, \code{\link{PACE}}.
#'
#' @export univBasisExpansion
#'
#' @examples
#' oldPar <- par(no.readonly = TRUE)
#' set.seed(1)
#'
#' # simulate some data (univariate)
#' sim <- simFunData(xVal = seq(-1,1, 0.01), M = 5, eFunType = "Poly",
#'                   eValType = "exponential", N = 7)$simData
#' noisy <- addError(sim, sd = 0.5) # a noisy version of sim
#'
#' simBasis <- univBasisExpansion(sim, pen = FALSE)
#' noiseBasis <- univBasisExpansion(noisy, pen = TRUE)
#'
#' # check reconstruction
#' par(mfrow = c(1,2))
#' plot(sim, main = "Simulated data")
#' plot(noisy, type = "p", pch = 20, cex = 0.3, add = TRUE)
#' legend("bottomright", c("sim", "noisy"), lty = c(1,NA), pch = c(NA, 20), pt.cex = 0.3)
#'
#' plot(funData(simBasis$functions@@xVal, simBasis$scores %*% simBasis$functions@@X),
#'      main = "Reconstruction based on\nbasis expansion")
#' plot(funData(noiseBasis$functions@@xVal, noiseBasis$scores %*% noiseBasis$functions@@X),
#'      lty = 3, add = TRUE)
#' legend("bottomright", c("sim, pen = FALSE", "noisy, pen = TRUE"), lty = c(1,3))
#'
#' par(oldPar)
univBasisExpansion <- function(funDataObject,
                               bs = "ps", # default: p-Splines
                               m = rep(2, dimSupp(funDataObject)), # gam/bam defaults
                               k = rep(-1, dimSupp(funDataObject)), # gam/bam defaults,
                               pen = FALSE) # use penalization to induce smoothing?
{
  if(dimSupp(funDataObject) > 2)
    stop("univBasisExpansion: Function is not implemented yet for objects of dimension > 2.")

  # special cases: k/m equal for all dimensions
  if(length(m) == 1)
    m <- rep(m, dimSupp(funDataObject))

  if(length(k) == 1)
    k <- rep(k, dimSupp(funDataObject))

  ### check correct dimensions
  if(any(c(length(m), length(k)) != dimSupp(funDataObject)))
    stop("univBasisExpansion: Specify 'm' and 'k' for each dimension of the functional data object!")

  if(length(bs) != 1)
    stop("univBasisExpansion: 'bs' must have length 1.")


  if(dimSupp(funDataObject) == 1)
  {
    if(pen)
      uniBasis <- splineBasis1Dpen(funDataObject, bs, m, k)
    else
      uniBasis <-  splineBasis1D(funDataObject, bs, m, k)
  }
  else
  {
    if(pen)
      uniBasis <- splineBasis2Dpen(funDataObject, bs, m, k)
    else
      uniBasis <-  splineBasis2D(funDataObject, bs, m, k)
  }

  return(uniBasis)
}


######## New stuff

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
#' principal component analysis)
#'
#' @param type A character string, specifying the basis for which the
#'   decomposition is to be calculated.
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
univExpansion <- function(type, scores, xVal, functions, params)
{
  params$scores <- scores

  if(is.numeric(xVal))
    xVal <- list(xVal)

  param$xVal <- xVal

  res <- switch(type,
                "uFPCA" = do.call(fpcaBasis, params),
                "splines1D" = do.call(splineFunction1D, params),
                "splines1Dpen" = do.call(splineFunction1D, params),
                "splines2D" = do.call(splineFunction2D, params),
                "splines2Dpen" = do.call(splineFunction2Dpen, params),
              #  "DCT2D" = ...,
                stop("Univariate Expansion for 'type' = ", type, " not defined!")
  )

  return(res$functions)
}

#' Calculate a linear combination of a functional principal component basis on
#' one-dimensional domains
#'
#' This function calculates  a linear combination of functional principal
#' component basis functions on one-dimensional domains.
#'
#' @param scores A matrix of dimension \eqn{N x K}, representing the \eqn{K}
#'   scores (coefficients) for each observation \eqn{i = 1, \ldots, N}.
#' @param xVal A list containing a vector of x-values.
#' @param functions A \code{funData} object, representing the FPCA basis.
#'
#' @return An object of class \code{funData} with \eqn{N} observations on
#'   \code{xVal}, corresponding to the linear combination of the functional
#'   principal components.
#'
#' @seealso univExpansion
fpcaFunction <- function(scores, xVal, functions)
{
  return(funData(xVal, scores %*% functions@X))
}


#' Calculate linear combinations of spline basis functions on one-dimensional
#' domains
#'
#' Given scores (coefficients), this function calculates a linear combination of
#' one-dimensional spline basis functions on one-dimensional domains based on the
#' \link[mgcv]{gam} function in the \pkg{mgcv} package.
#'
#' @param scores A matrix of dimension \eqn{N x K}, representing the \eqn{K}
#'   scores (coefficients) for each observation \eqn{i = 1, \ldots, N}.
#' @param xVal A list containing a vector of x-values.
#' @param bs A character string, specifying the type of basis functions to be
#'   used. Please refer to \code{\link[mgcv]{smooth.terms}} for a list of
#'   possible basis functions.
#' @param m A numeric, the order of the spline basis. See  \code{\link[mgcv]{s}}
#'   for details.
#' @param k A numeric, the number of basis functions used.  See
#'   \code{\link[mgcv]{s}} for details.
#'
#' @return An object of class \code{funData} with \eqn{N} observations on \code{xVal},
#'   corresponding to the linear combination of spline basis functions.
#'
#'  @seealso univExpansion
#'
#' @importFrom mgcv gam
splineFunction1D <- function(scores, xVal, bs, m, k)
{
  N <- nrow(scores)

  x <- xVal[[1]]

  # spline design matrix via gam
  desMat <- mgcv::gam(rep(0, length(x)) ~ s(x, bs = bs, m = m, k = k), fit = FALSE)$X

  # calculate functions as linear combination of splines
  res <- funData(xVal,
                 scores %*% t(desMat))

  return(res)
}


#' Calculate linear combinations of spline basis functions on two-dimensional
#' domains
#'
#' Given scores (coefficients), this function calculates a linear combination of
#' two-dimensional spline tensor basis functions on one-dimensional domains
#' based on the \link[mgcv]{gam} function in the \pkg{mgcv} package. If the
#' scores have been calculated based on a penalized tensor spline basis, use
#' \code{splineFunction2Dpen} instead (which runs \link[mgcv]{bam} instead of
#' \link[mgcv]{gam}).
#'
#' @param scores A matrix of dimension \eqn{N x K}, representing the \eqn{K}
#'   scores (coefficients) for each observation \eqn{i = 1, \ldots, N}.
#' @param xVal A list containing a two numeric vectors, corresponding to the x-
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
#'   tow-dimensional domain specified by \code{xVal}, corresponding to the
#'   linear combination of spline basis functions.
#'
#'  @seealso univExpansion
#'
#' @importFrom mgcv gam
splineFunction2D <- function(scores, xVal, bs, m, k)
{
  N <- nrow(scores)

  coord <- expand.grid(x = xVal[[1]], y = xVal[[2]])

  # spline design matrix via gam
  desMat <- mgcv::gam(rep(0, dim(coord)[1]) ~ te(coord$x, coord$y, bs = bs, m = m, k = k), data = coord, fit = FALSE)$X

  # calculate functions as linear combination of splines
  res <- funData(xVal,
                 array(scores %*% t(desMat),
                       dim = c(N, length(xVal[[1]]), length(xVal[[2]]))))

  return(res)
}


#' @rdname splineFunction2D
#'
#' @importFrom mgcv gam
splineFunction2Dpen <- function(scores, xVal, bs, m, k)
{
  N <- nrow(scores)

  coord <- expand.grid(x = xVal[[1]], y = xVal[[2]])

  # spline design matrix via gam
  desMat <- mgcv::bam(rep(0, dim(coord)[1]) ~ te(coord$x, coord$y, bs = bs, m = m, k = k), data = coord, fit = FALSE)$X

  # calculate functions as linear combination of splines
  res <- funData(xVal,
                 array(scores %*% t(desMat),
                       dim = c(N, length(xVal[[1]]), length(xVal[[2]]))))

  return(res)
}
