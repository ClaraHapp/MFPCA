## functional PCA: consider X-values as additional input
# slightly adapted version of fpca.sc-function in refund-package
# Internal function, called by PACE in context of funData objects

#' Calculate univariate functional PCA
#' 
#' This function is a slightly adapted version of the 
#' \code{\link[refund]{fpca.sc}} function in the \pkg{refund} package for 
#' calculating univariate functional principal components based on a smoothed 
#' covariance function. The smoothing basis functions are penalized splines.
#' 
#' @param X A vector of xValues.
#' @param Y A matrix of observed functions (by row).
#' @param Y.pred A matrix of functions (by row) to be approximated using the 
#'   functional principal components. Defaults to \code{NULL}, i.e. the 
#'   prediction is made for the functions in \code{Y}.
#' @param nbasis An integer, giving the number of B-spline basis to use. 
#'   Defaults to \code{10}.
#' @param pve A value between 0 and 1, giving the percentage of variance 
#'   explained in the data by the functional principal components. This value is
#'   used to choose the number of principal components. Defaults to \code{0.99}
#' @param npc The number of principal components to be estimated. Defaults to 
#'   \code{NULL}. If given, this overrides \code{pve}.
#' @param makePD Logical, should positive definiteness be enforced for the 
#'   covariance estimate? Defaults to \code{FALSE}.
#' @param cov.weight.type The type of weighting used for the smooth covariance
#'   estimate. Defaults to \code{"none"}, i.e. no weighting. Alternatively, 
#'   \code{"counts"} (corresponds to \code{\link[refund]{fpca.sc}} ) weights the pointwise estimates of the covariance function
#'   by the number of observation points.
#'   
#' @return \item{fit}{The approximation of \code{Y.pred} (if \code{NULL}, the 
#'   approximation of \code{Y}) based on the functional principal components.} 
#'   \item{scores}{A matrix containing the estimated scores (observations by 
#'   row).} \item{mu}{The estimated mean function.} \item{efunctions}{A matrix 
#'   containing the estimated eigenfunctions (by row).} \item{evalues}{The 
#'   estimated eigenvalues.} \item{npc}{The number of principal comopnents that 
#'   were calculated.} \item{sigma2}{The estimated variance of the measurement 
#'   error.}  \item{estVar}{The estimated smooth variance function of the data.}
#'   
#' @seealso \code{\link[refund]{fpca.sc}}, \code{\link{PACE}}
#'   
#' @references Di, C., Crainiceanu, C., Caffo, B., and Punjabi, N. (2009). 
#'   Multilevel functional principal component analysis. Annals of Applied 
#'   Statistics, 3, 458--488. Yao, F., Mueller, H.-G., and Wang, J.-L. (2005). 
#'   Functional data analysis for sparse longitudinal data. Journal of the 
#'   American Statistical Association, 100, 577--590.
#'   
#' @importFrom mgcv gam predict.gam s te
#'   
#' @keywords internal
.PACE <- function(X, Y, Y.pred = NULL, nbasis = 10, pve = 0.99, npc = NULL, makePD = FALSE, cov.weight.type = "none")
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
  cov.weights <- switch(cov.weight.type,
                        none = rep(1, D^2),
                        counts = as.vector(cov.count),
                        stop("cov.weight.type ", cov.weight.type, " unknown in smooth covariance estimation"))
  
  npc.0 = matrix(mgcv::predict.gam(mgcv::gam(as.vector(G.0)~te(row.vec, col.vec, k = nbasis), weights = cov.weights),
                                   newdata = data.frame(row.vec = row.vec, col.vec = col.vec)), D, D)
  npc.0 = (npc.0 + t(npc.0))/2
  # no extra-option (useSymm) as in fpca.sc-method
  if (makePD) { # see fpca.sc
    npc.0 <- {
      tmp <- Matrix::nearPD(npc.0, corr = FALSE, keepDiag = FALSE,
                            do2eigen = TRUE, trace = options()$verbose)
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
  sigma2 <- max(1/(X[T1.max]-X[T1.min]) * sum(DIAG*w, na.rm = TRUE), 0) #max(1/T.len * sum(DIAG*w), 0)
  ####
  D.inv = diag(1/evalues, nrow = npc, ncol = npc)
  Z = efunctions
  Y.tilde = Y.pred - matrix(mu, I.pred, D, byrow = TRUE)
  fit = matrix(0, nrow = I.pred, ncol = D)
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
    scores[i.subj, ] = ZtZ_sD.inv %*% crossprod(Zcur, Y.tilde[i.subj, obs.points])
    fit[i.subj, ] = t(as.matrix(mu)) + tcrossprod(scores[i.subj, ], efunctions)
  }
  ret.objects = c("fit", "scores", "mu", "efunctions", "evalues",
                  "npc", "sigma2") # add sigma2 to output
  ret = lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret) = ret.objects
  ret$estVar <- diag(cov.hat)
  return(ret)
}

#' Univariate functional principal component analysis by smoothed covariance
#' 
#' This function calculates a univariate functional principal components 
#' analysis by smoothed covariance based on code from 
#' \code{\link[refund]{fpca.sc}} (package \pkg{refund}).
#' 
#' @section Warning: This function works only for univariate functional data 
#'   observed on one-dimensional domains.
#'   
#' @param funDataObject An object of class \code{\link[funData]{funData}} or 
#'   \code{\link[funData]{irregFunData}} containing the functional data 
#'   observed, for which the functional principal component analysis is 
#'   calculated. If the data is sampled irregularly (i.e. of class 
#'   \code{\link[funData]{irregFunData}}), \code{funDataObject} is transformed 
#'   to a \code{\link[funData]{funData}} object first.
#' @param predData  An object of class \code{\link[funData]{funData}}, for which
#'   estimated trajectories based on a truncated Karhunen-Loeve representation 
#'   should be estimated. Defaults to \code{NULL}, which implies prediction for 
#'   the given data.
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
#' @param cov.weight.type The type of weighting used for the smooth covariance 
#'   estimate. Defaults to \code{"none"}, i.e. no weighting. Alternatively, 
#'   \code{"counts"} (corresponds to \code{\link[refund]{fpca.sc}} ) weights the
#'   pointwise estimates of the covariance function by the number of observation
#'   points.
#'   
#' @return \item{mu}{A \code{\link[funData]{funData}} object with one 
#'   observation, corresponding to the mean function.} \item{values}{A vector 
#'   containing the estimated eigenvalues.} \item{functions}{A 
#'   \code{\link[funData]{funData}} object containing the estimated functional 
#'   principal components.} \item{scores}{An matrix of estimated scores for the 
#'   observations in \code{funDataObject}. Each row corresponds to the scores of
#'   one observation.} \item{fit}{A \code{\link[funData]{funData}} object 
#'   containing the estimated trajectories based on the truncated Karhunen-Loeve
#'   representation and the estimated scores and functional principal components
#'   for \code{predData} (if this is not \code{NULL}) or \code{funDataObject} 
#'   (if \code{predData} is \code{NULL}).} \item{npc}{The number of functional 
#'   principal components: either the supplied \code{npc}, or the minimum number
#'   of basis functions needed to explain proportion \code{pve} of the variance 
#'   in the observed curves (cf. \code{\link[refund]{fpca.sc}}).} 
#'   \item{sigma2}{The estimated measurement error variance (cf. 
#'   \code{\link[refund]{fpca.sc}}).} \item{estVar}{The estimated smooth
#'   variance function of the data.}
#'   
#' @seealso \code{\link[funData]{funData}}, \code{\link[refund]{fpca.sc}}, 
#'   \code{\link{fpcaBasis}}, \code{\link{univDecomp}}
#'   
#' @export PACE
#'   
#' @examples
#' \donttest{
#'   oldPar <- par(no.readonly = TRUE)
#' 
#'   # simulate data
#'   sim <- simFunData(argvals = seq(-1,1,0.01), M = 5, eFunType = "Poly",
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
PACE <- function(funDataObject, predData = NULL, nbasis = 10, pve = 0.99, npc = NULL, makePD = FALSE, cov.weight.type = "none")
{
  if(dimSupp(funDataObject) != 1)
    stop("PACE: Implemented only for funData objects with one-dimensional support.")
  
  if(methods::is(funDataObject, "irregFunData")) # for irregular functional data, use funData representation
    funDataObject <- as.funData(funDataObject)
  
  if(!is.null(predData))
  {
    if(!isTRUE(all.equal(funDataObject@argvals, predData@argvals)))
      stop("PACE: funDataObject and predData must be defined on the same domains!")
    
    Y.pred = predData@X
  }
  else
  {
    Y.pred = NULL # use only funDataObject
  }
  
  res <- .PACE(X = funDataObject@argvals[[1]], funDataObject@X, Y.pred = Y.pred,
               nbasis = nbasis, pve = pve, npc = npc, makePD = makePD,
               cov.weight.type = cov.weight.type)
  
  return(list(mu = funData(funDataObject@argvals, matrix(res$mu, nrow = 1)),
              values = res$evalues,
              functions = funData(funDataObject@argvals, t(res$efunctions)),
              scores = res$scores,
              fit = funData(funDataObject@argvals, res$fit),
              npc = res$npc,
              sigma2 = res$sigma2,
              estVar = funData(funDataObject@argvals, matrix(res$estVar, nrow = 1))
  ))
}
