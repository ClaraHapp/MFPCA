# define global variable i, used by the foreach package and confusing R CMD CHECK
globalVariables('i')

#' @import funData
NULL

#' Calculate an unpenalized spline basis representation for functional data on
#' one-dimensional domains
#'
#' This function calculates an unpenalized spline basis representation for
#' functional data on one-dimensional domains based on the \link[mgcv]{gam}
#' function in the \pkg{mgcv} package.
#'
#' @param funDataObject An object of class \code{\link[funData]{funData}}
#'   containing the observed functional data samples and for which the bases
#'   representation is calculated.
#' @param bs The type of basis functions to be used. Please refer to
#'   \code{link[mgcv]{smooth.terms}} for a list of possible basis functions.
#' @param m The order of the spline basis. See  \code{link[mgcv]{s}} for
#'   details.
#' @param k The number of basis functions used.  See  \code{link[mgcv]{s}} for
#'   details.
#'
#' @return scores A matrix of weights with dimension \code{N x k}, reflecting
#'   the weights for each basis function in each observation.
#' @return function A \code{\link[funData]{funData}} object representing the
#'   basis functions.
#'
#'   @importFrom mgcv gam
#'
#' @keywords internal
splineBasis1D <- function(funDataObject, bs, m, k)
{
  N <- nObs(funDataObject)

  x <- funDataObject@xVal[[1]]

  # spline design matrix via gam
  desMat <- mgcv::gam(funDataObject@X[1, ] ~ s(x, bs = bs, m = m, k = k), fit = FALSE)$X

  # weights via lm -> no penalization
  scores <- t(apply(funDataObject@X, 1, function(f, dM){lm(f ~ dM - 1)$coef}, dM = desMat))

  return(list(scores = scores,
              functions = funData(x, t(desMat))))
}

#' Calculate an unpenalized spline basis representation for functional data on
#' two-dimensional domains
#'
#' This function calculates an unpenalized tensor product spline basis
#' representation for functional data on two-dimensional domains based on the
#' \link[mgcv]{gam} function in the \pkg{mgcv} package.
#'
#' @param funDataObject An object of class \code{\link[funData]{funData}}
#'   containing the observed functional data samples and for which the bases
#'   representation is calculated.
#' @param bs An array (or a single character), the type of basis functions to be
#'   used. Please refer to \code{link[mgcv]{te}} for a list of possible basis
#'   functions.
#' @param m An array (or a single character), the order of the spline basis. See
#'   \code{link[mgcv]{s}} for details.
#' @param k An array (or a single character), the number of basis functions
#'   used.  See  \code{link[mgcv]{s}} for details.
#'
#' @return scores A matrix of weights with dimension \code{N x k}, reflecting
#'   the weights for each basis function in each observation.
#' @return function A \code{\link[funData]{funData}} object representing the
#'   basis functions.
#' @return basisLong A matrix representing the basis functions in vectorized
#'   form.
#'
#'@importFrom foreach %do%
#'   @importFrom mgcv gam
#'
#' @keywords internal
splineBasis2D <- function(funDataObject, bs, m, k)
{
  N <- nObs(funDataObject)

  x <- rep(funDataObject@xVal[[1]], times = length(funDataObject@xVal[[2]]))
  y <- rep(funDataObject@xVal[[2]], each = length(funDataObject@xVal[[1]]))

  # spline design matrix via gam
  desMat <- mgcv::gam(as.vector(funDataObject@X[1,,]) ~ te(x, y, bs = bs, m = m, k = k), fit = FALSE)$X

  # weights via lm -> no penalization
  scores <- t(apply(funDataObject@X, 1, function(f, dM){lm(as.vector(f) ~ dM - 1)$coef}, dM = desMat))

  # extract basis functions (in the correct dimensions)
  tmp <-foreach::foreach(i = 1:dim(desMat)[2], .combine = function(x, y){abind::abind(x, y, along = 3)})%do%{
    matrix(desMat[, i], nrow = length(funDataObject@xVal[[1]]), ncol = length(funDataObject@xVal[[2]]))}

  return(list(scores = scores,
              functions = funData(funDataObject@xVal, aperm(tmp, c(3,1,2)) ),
              basisLong = desMat))
}




#' Calculate a penalized spline basis representation for functional data on
#' one-dimensional domains
#'
#' This function calculates a penalized spline basis representation for
#' functional data on one-dimensional domains based on the \link[mgcv]{gam}
#' function in the \pkg{mgcv} package.
#'
#' @param funDataObject An object of class \code{\link[funData]{funData}}
#'   containing the observed functional data samples and for which the bases
#'   representation is calculated.
#' @param bs The type of basis functions to be used. Please refer to
#'   \code{link[mgcv]{smooth.terms}} for a list of possible basis functions.
#' @param m The order of the spline basis. See  \code{link[mgcv]{s}} for
#'   details.
#' @param k The number of basis functions used.  See  \code{link[mgcv]{s}} for
#'   details.
#'
#' @return scores A matrix of weights with dimension \code{N x k}, reflecting
#'   the weights for each basis function in each observation.
#' @return function A \code{\link[funData]{funData}} object representing the
#'   basis functions.
#'
#' @importFrom foreach %do%
#' @importFrom mgcv gam
#'
#' @keywords internal
splineBasis1Dpen <- function(funDataObject, bs, m, k)
{
  N <- nObs(funDataObject)

  x <- funDataObject@xVal[[1]]

  scores <- foreach::foreach(i = 1:(N-1), .combine = "rbind")%do%{
    g <- mgcv::gam(funDataObject@X[i, ] ~ s(x, bs = bs, m = m, k = k), method = "REML")
    g$coef
  }

  # last extra to extract model matrix
  g <- mgcv::gam(funDataObject@X[N, ] ~ s(x, bs = bs, m = m, k = k), method = "REML")

  scores <- rbind(scores, g$coef)

  return(list(scores = scores,
              functions = funData(x, t(model.matrix(g)))))
}


#' Calculate a penalized spline basis representation for functional data on
#' two-dimensional domains
#'
#' This function calculates a penalized tensor product spline basis
#' representation for functional data on two-dimensional domains based on the
#' \link[mgcv]{bam} function in the \pkg{mgcv} package (for large GAMs).
#'
#' @param funDataObject An object of class \code{\link[funData]{funData}}
#'   containing the observed functional data samples and for which the bases
#'   representation is calculated.
#' @param bs An array (or a single character), the type of basis functions to be
#'   used. Please refer to \code{link[mgcv]{te}} for a list of possible basis
#'   functions.
#' @param m An array (or a single character), the order of the spline basis. See
#'   \code{link[mgcv]{s}} for details.
#' @param k An array (or a single character), the number of basis functions
#'   used.  See  \code{link[mgcv]{s}} for details.
#'
#' @return scores A matrix of weights with dimension \code{N x k}, reflecting
#'   the weights for each basis function in each observation.
#' @return function A \code{\link[funData]{funData}} object representing the
#'   basis functions.
#' @return basisLong A matrix representing the basis functions in vectorized
#'   form.
#'
#'@importFrom foreach %do%
#'   @importFrom mgcv bam
#'
#' @keywords internal
splineBasis2Dpen <- function(funDataObject, bs, m, k)
{
  N <- nObs(funDataObject)

  x <- rep(funDataObject@xVal[[1]], times = length(funDataObject@xVal[[2]]))
  y <- rep(funDataObject@xVal[[2]], each = length(funDataObject@xVal[[1]]))

  scores <- foreach::foreach(i = 1:(N-1), .combine = "rbind")%do%{
    g <- mgcv::bam(as.vector(funDataObject@X[i, , ]) ~ te(x, y, bs = bs, m = m, k = k), method = "REML")
    g$coef
  }

  # fit the last one extra in order to extract model matrix
  g <- mgcv::bam(as.vector(funDataObject@X[N, , ]) ~ te(x, y, bs = bs, m = m, k = k), method = "REML")

  scores <- rbind(scores, g$coef)

  basisLong <- model.matrix(g)

  # extract basis functions (in the correct dimensions)
  tmp <-foreach::foreach(i = 1:prod(k), .combine = function(x, y){abind::abind(x, y, along = 3)})%do%{
    matrix(basisLong[, i], nrow = length(funDataObject@xVal[[1]]), ncol = length(funDataObject@xVal[[2]]))}

  return(list(scores = scores,
              functions = funData(funDataObject@xVal, aperm(tmp, c(3,1,2))),
              basisLong = basisLong))
}

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
univExpansion(type, scores, functions, params)
{
  res <- switch(type,
                "uFPCA" = ...,
                "splines1D" = ...,
                "splines1Dpen" = ...,
                "splines2D" = ...,
                "splines2Dpen" = ...,
                "DCT2D" = ...,
                stop("Univariate Decomposition for 'type' = ", type, " not defined!")
  )

  return(res$functions)
}


