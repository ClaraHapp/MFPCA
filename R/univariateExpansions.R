# define global variable i, used by the foreach package and confusing R CMD CHECK
globalVariables('i')

#' @import funData
NULL

## functional PCA: consider X-values as additional input
# slightly adapted version of fpca.sc-function in refund-package
# Internal function, called by PACE in context of funData objects

#' @importFrom mgcv gam predict.gam
.PACE <- function(X, Y, Y.pred = NULL, nbasis = 10, pve = 0.99, npc = NULL, makePD = FALSE)
{
  if (is.null(Y.pred))
    Y.pred = Y
  D = NCOL(Y)
  if(D != length(X)) # check if number of observation points in X & Y are identical
    stop("different number of (potential) observation points differs in X and Y!")
  I = NROW(Y)
  I.pred = NROW(Y.pred)
  d.vec = rep(X, each = I) # use given X-values for estimation of mu
  gam0 = mgcv::gam(as.vector(Y) ~ s(d.vec, k = nbasis))
  mu = mgcv::predict.gam(gam0, newdata = data.frame(d.vec = X))
  Y.tilde = Y - matrix(mu, I, D, byrow = TRUE)
  cov.sum = cov.count = cov.mean = matrix(0, D, D)
  for (i in 1:I) {
    obs.points = which(!is.na(Y[i, ]))
    cov.count[obs.points, obs.points] = cov.count[obs.points, obs.points] + 1
    cov.sum[obs.points, obs.points] = cov.sum[obs.points, obs.points] + tcrossprod(Y.tilde[i, obs.points])
  }
  G.0 = ifelse(cov.count == 0, NA, cov.sum/cov.count)
  diag.G0 = diag(G.0)
  diag(G.0) = NA
  row.vec = rep(X, each = D) # use given X-values
  col.vec = rep(X, D) # use given X-values
  npc.0 = matrix(mgcv::predict.gam(mgcv::gam(as.vector(G.0)~te(row.vec, col.vec, k = nbasis)),
                         newdata = data.frame(row.vec = row.vec, col.vec = col.vec)), D, D)
  npc.0 = (npc.0 + t(npc.0))/2
  # no extra-option (useSymm) as in fpca.sc-method
  if (makePD) { # see fpca.sc
    npc.0 <- {
      tmp <- Matrix::nearPD(npc.0, corr = FALSE, keepDiag = FALSE,
                            do2eigen = TRUE, trace = TRUE)
      as.matrix(tmp$mat)
    }
  }

  ### numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Chapter 8)
  w <- funData:::.intWeights(X, method = "trapezoidal")
  Wsqrt <- diag(sqrt(w))
  Winvsqrt <- diag(1/(sqrt(w)))
  V <- Wsqrt %*% npc.0 %*% Wsqrt
  evalues = eigen(V, symmetric = TRUE, only.values = TRUE)$values
  ###
  evalues = replace(evalues, which(evalues <= 0), 0)
  npc = ifelse(is.null(npc), min(which(cumsum(evalues)/sum(evalues) > pve)), npc)
  efunctions = matrix(Winvsqrt%*%eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = D, ncol = npc)
  evalues = eigen(V, symmetric = TRUE, only.values = TRUE)$values[1:npc]  # use correct matrix for eigenvalue problem
  cov.hat = efunctions %*% tcrossprod(diag(evalues, nrow = npc, ncol = npc), efunctions)
  ### numerical integration for estimation of sigma2
  T.len <- X[D] - X[1] # total interval length
  T1.min <- min(which(X >= X[1] + 0.25*T.len)) # left bound of narrower interval T1
  T1.max <- max(which(X <= X[D] - 0.25*T.len)) # right bound of narrower interval T1
  DIAG = (diag.G0 - diag(cov.hat))[T1.min :T1.max] # function values
  # weights
  w <- funData:::.intWeights(X[T1.min:T1.max], method = "trapezoidal")
  sigma2 <- max(1/(X[T1.max]-X[T1.min]) * sum(DIAG*w), 0) #max(1/T.len * sum(DIAG*w), 0)
  ####
  D.inv = diag(1/evalues, nrow = npc, ncol = npc)
  Z = efunctions
  Y.tilde = Y.pred - matrix(mu, I.pred, D, byrow = TRUE)
  Yhat = matrix(0, nrow = I.pred, ncol = D)
  scores = matrix(NA, nrow = I.pred, ncol = npc)
  # no calculation of confidence bands, no variance matrix
  for (i.subj in 1:I.pred) {
    obs.points = which(!is.na(Y.pred[i.subj, ]))
    if (sigma2 == 0 & length(obs.points) < npc) {
      stop("Measurement error estimated to be zero and there are fewer observed points than PCs; scores cannot be estimated.")
    }
    Zcur = matrix(Z[obs.points, ], nrow = length(obs.points),
                  ncol = dim(Z)[2])
    ZtZ_sD.inv = solve(crossprod(Zcur) + sigma2 * D.inv)
    scores[i.subj, ] = ZtZ_sD.inv %*% t(Zcur) %*% (Y.tilde[i.subj, obs.points])
    Yhat[i.subj, ] = t(as.matrix(mu)) + scores[i.subj, ] %*%
      t(efunctions)
  }
  ret.objects = c("Yhat", "scores", "mu", "efunctions", "evalues",
                  "npc", "sigma2") # add sigma2 to output
  ret = lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret) = ret.objects
  return(ret)
}

#' Univariate functional principal component analysis by smoothed covariance
#'
#' This function calculates a univariate functional principal components
#' analysis by smoothed covariance based on code from \link[refund]{fpca.sc}
#' (pacakge \pkg{refund}).
#'
#' @section Warning: This function works only for univariate functional data
#'   observed on one-dimensional domains.
#'
#' @param funDataObject An object of class \code{\link[funData]{funData}}
#'   containing the functional data observed, for which the functional principal
#'   component analysis is calculated.
#' @param predData  An object of class \code{\link[funData]{funData}}, for which
#'   estimated trajectories based on a truncated Karhunen-Lo\`{e}ve representation
#'   should be estimated. Defaults to \code{NULL}, which implies prediction for
#'   the given data.
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
#' @return \item{fit}{A \code{\link[funData]{funData}} object containing the
#'   estimated trajectories based on the truncated Karhunen-Lo\`{e}ve representation
#'   and the estimated scores and functional principal components for
#'   \code{predData} (if this is not \code{NULL}) or \code{funDataObject} (if
#'   \code{predData} is \code{NULL}).} \item{scores}{An matrix of estimated
#'   scores for the observations in \code{funDataObject}. Each row corresponds
#'   to the scores of one observation.} \item{mu}{A
#'   \code{\link[funData]{funData}} object with one observation, corresponding
#'   to the mean function.} \item{functions}{A \code{\link[funData]{funData}}
#'   object containing the estimated functional principal components.}
#'   \item{values}{A vector containing the estimated eigenvalues.}
#'   \item{npc}{The number of functional principal components: either the
#'   supplied \code{npc}, or the minimum number of basis functions needed to explain
#'   proportion \code{pve} of the variance in the observed curves (cf.
#'   \code{\link[refund]{fpca.sc}}).} \item{sigma2}{The estimated measurement
#'   error variance (cf. \code{\link[refund]{fpca.sc}}).}
#'
#' @seealso \code{\link[funData]{funData}}, \code{\link[refund]{fpca.sc}}.
#'
#' @export PACE
#'
#' @examples
#' \donttest{
#'   oldPar <- par(no.readonly = TRUE)
#'
#'   # simulate data
#'   sim <- simFunData(xVal = seq(-1,1, 0.01), M = 5, eFunType = "Poly",
#'                     eValType = "exponential", N = 100)
#'
#'   # calculate univariate FPCA
#'   pca <- PACE(sim$simData, npc = 5)
#'
#'   # Plot the results
#'   par(mfrow = c(1,2))
#'   plot(sim$trueFuns, lwd = 2, main = "Eigenfunctions")
#'   # flip estimated functions for correct signs
#'   plot(flipFuns(sim$trueFuns,pca$functions), lty = 2, add = TRUE)
#'   legend("bottomright", c("True", "Estimate"), lwd = c(2,1), lty = c(1,2))
#'
#'   plot(sim$simData, lwd = 2, main = "Some Observations", obs = 1:7)
#'   plot(pca$fit, lty = 2, obs = 1:7, add = TRUE) # estimates are almost equal to true values
#'   legend("bottomright", c("True", "Estimate"), lwd = c(2,1), lty = c(1,2))
#'
#'   par(oldPar)
#' }
PACE <- function(funDataObject, predData = NULL, nbasis = 10, pve = 0.99, npc = NULL, makePD = FALSE)
{
  if(dimSupp(funDataObject) != 1)
    stop("PACE: Implemented only for funData objects with one-dimensional support.")

  if(!is.null(predData))
  {
    if(!isTRUE(all.equal(funDataObject@xVal, predData@xVal)))
      stop("PACE: funDataObject and predData must be defined on the same domains!")

    Y.pred = predData@X
  }
  else
  {
    Y.pred = NULL # use only funDataObject
  }

  res <- .PACE(X = funDataObject@xVal[[1]], funDataObject@X, Y.pred = Y.pred, nbasis = nbasis, pve = pve, npc = npc, makePD = makePD)

  return(list(fit = funData(funDataObject@xVal, res$Yhat),
              scores = res$scores,
              mu = funData(funDataObject@xVal, matrix(res$mu, nrow = 1)),
              functions = funData(funDataObject@xVal, t(res$efunctions)),
              values = res$evalues,
              npc = res$npc,
              sigma2 = res$sigma2
  ))
}


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
#'   (\code{"ps"}) for functions on one-dimensional domain and to thin plate
#'   regression splines (\code{"tp"}) in case of two-dimensional support. Please
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
                               bs = NULL, # default (see later) 1D: p-Splines, 2D: Tensor Product Splines (Thin Plate)
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

  if(is.null(bs)) # default
    bs <- ifelse(dimSupp(funDataObject)==1, "ps", "tp")

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
