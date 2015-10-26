#' @import funData
NULL

#' Utility function that calculates matrix of basis-scalar products (one dimension)
#'
#' If the element \eqn{X^{(j)}}{X^{(j)}} is expanded in basis functions \eqn{b_i^{(j)}(t),~ i = 1, \ldots, K_j}{b_i(t)},
#' this function calculates the \eqn{K_j \times K_j}{K_j  \times K_j} matrix \eqn{B^{(jj)}}{B^{(jj)}} with entries
#' \deqn{B^{(jj)}_{mn} = \int_{\mathcal{T_j}} b_m^{(j)}(t) b_n^{(j)}(t) \mathrm{d} t}.
#'
#' @section Warning: This function is implemented only for functions on one- or two-dimensional domains.
#'
#' @param basisFunctions Array of \code{npc} basis functions of dimensions \code{npc x M1} or \code{npc x M1 x M2}.
#' @param dimSupp dimension of the support of the basis functions (1 or 2)
#' @param xVal List of corresponding x-values.
#'
#' @return A matrix containing the scalar product of all combinations of basis functions (matrix \eqn{B^{(j)}})
#'
#' @seealso \code{\link{MFPCA}}, \code{\link[funData]{dimSupp}}.
#'
#' @keywords internal
.calcBasisIntegrals <- function(basisFunctions, dimSupp, xVal)
{
  npc <- dim(basisFunctions)[1]

  #  integral basis matrix
  B <- array(0, dim = c(npc, npc))


  if(dimSupp == 1) # one-dimensional domain
  {
    w <- funData:::.intWeights(xVal[[1]])

    for(m in 1:npc)
    {
      for(n in 1:m)
        B[m, n] <- B[n, m] <- (basisFunctions[m, ]* basisFunctions[n, ])%*%w
    }
  }
  else # two-dimesional domain (otherwise function stops before!)
  {
    w1 <- t(funData:::.intWeights(xVal[[1]]))
    w2 <- funData:::.intWeights(xVal[[2]])

    for(m in 1:npc)
    {
      for(n in 1:m)
        B[m, n] <- B[n, m] <-  w1 %*%(basisFunctions[m, , ]* basisFunctions[n, ,])%*%w2
    }
  }

  return(B)
}


#' Multivariate Fuctional Principal Component Analysis for Functions on
#' Different (Dimensional) Domains
#'
#' This function calculates a multivariate fuctional principal component
#' analysis (MFPCA) based on iid. observations \eqn{x_1, \ldots, x_N} of a
#' multivariate functional data-generating process \eqn{X = (X^{(1)}, \ldots
#' X^{(p)})}{X = X^(1), \ldots, X^(p)} with elements \eqn{X^{(j)} \in
#' L^2(\mathcal{T}_j)}{X^(j) in L^2(calT_j)} defined on a domain
#' \eqn{\mathcal{T}_j \subset IR^{d_j}}{calT_j of IR^{d_j}}. In particular, the
#' elements can be defined on different (dimensional) domains. The results
#' contain the estimated multivariate functional principal components \eqn{\hat
#' \psi_1, \ldots, \hat \psi_M} (having the same structure as \eqn{x_i}), the
#' associated eigenvalues \eqn{\hat \nu_1 \geq \ldots \geq \hat \nu_M > 0} and
#' the individual scores \eqn{\hat \rho_{im} = \widehat{<x_i, \psi_m>}}{\hat
#' \rho_{im} = \hat{<x_i, \psi_m>}}. Moreover, estimated trajectories for each
#' observation based on the truncated Karhunen-Lo\`{e}ve representation
#' \deqn{\hat x_i = \sum_{m = 1}^M \hat \rho_{im} \hat \psi_m}{\hat x_i =
#' \sum_{m = 1}^M \hat \rho_{im} \hat \psi_m} are given. The implementation of
#' the  observations \eqn{x_i = (x_i^{(1)}, \ldots , x_i^{(p)}),~ i = 1 ,
#' \ldots, N}{x_i = (x_i^(1), \ldots , x_i^(p)), i = 1 , \ldots, N} and
#' multivariate functional principal components \eqn{\hat \psi_1, \ldots, \hat
#' \psi_M} uses the \code{\link[funData]{multiFunData}} class, which is defined
#' in the package \pkg{funData}.
#'
#' \subsection{Weighted MFPCA:}{If the elements vary considerably in domain,
#' range or variation, a weight vector \eqn{w_1 , \ldots, w_p} can be supplied
#' and the MFPCA is based on the weighted scalar product \deqn{<<f,g>>_w =
#' \sum_{j = 1}^p w_j \int_{\mathcal{T_j}} f^{(j)}(t) g^{(j)}(t) \mathrm{d}
#' t}{<<f,g>>_w = \sum_{j = 1}^p w_j \int_\calT_j f^(j)(t) g^(j)(t) d t} and the
#' corresponding weighted covariance operator \eqn{\Gamma_w}.}
#'
#' \subsection{Pointwise bootstrap confidence bands:}{Optionally, 95\% pointwise
#' bootstrap confidence bands are generated for the multivariate functional
#' principal components \eqn{\hat \psi_1, \ldots, \hat \psi_M}{\hat \psi_1,
#' \ldots, \hat \psi_M}.}
#'
#'
#' \subsection{Univariate Expansions:}{The multivariate functional principal
#' component analysis relies on a univariate basis expansion for each element
#' \eqn{X^{(j)}}{X^(j)}. It can be supplied in several forms: \itemize{ \item
#' Univariate functional principal component analysis. Then \code{uniExpansions
#' = list(type = "uFPCA", nbasis, pve, npc, makePD)}, where \code{nbasis, pve,
#' npc, makePD} are parameters passed to the \code{\link{PACE}} function for
#' calculating the univariate functional principal component analysis. \item
#' Spline basis functions (not penalized). Then \code{uniExpansions = list(type
#' = "splines", bs, m, k)}, where \code{bs, m, k} are passed to the function
#' \code{\link{univBasisExpansion}}. \item Spline basis functions (with
#' smoothness penalty). Then \code{uniExpansions = list(type = "splinesPen", bs,
#' m, k)}, where \code{bs, m, k} are passed to the function
#' \code{\link{univBasisExpansion}}. \item General basis functions. Then
#' \code{uniExpansions = list(functions, scores)}, where \code{functions} is a
#' \code{\link[funData]{multiFunData}} object containing the basis functions and
#' \code{scores} is an array of dimensions \code{N x M_1 x \ldots x M_j} with
#' the corresponding basis weights.}}
#'
#' @param mFData A  \code{\link[funData]{multiFunData}} object containing the
#'   observations \eqn{x_i = (x_i^{(1)}, \ldots , x_i^{(p)}),~ i = 1 , \ldots,
#'   N}{x_i = (x_i^(1), \ldots , x_i^(p)), i = 1 , \ldots, N}.
#' @param M The number of multivariate functional principal components to
#'   calculate.
#' @param uniExpansions A list characterizing the (univariate) expansion that is
#'   calculated for each element. See Details.
#' @param weights An optional vector of weights, defaults to 1 for each element.
#'   See Details.
#' @param Yhat Logical. If \code{TRUE} a truncated multivariate
#'   Karhunen-Lo\`{e}ve represenation for the data is calcualted based on the
#'   estimated scores and eigenfunctions.
#' @param bootstrap Logical. If \code{TRUE}, pointwise bootstrap confidence
#'   bands are calculated for the multivariate functional principal components.
#'   Defaults to \code{FALSE}. See Details.
#' @param nBootstrap The number of bootstrap iterations to use. Defaults to
#'   \code{NULL}, which leads to an error, if \code{bootstrap = TRUE}.
#' @param bootstrapAlpha A vector of numerics (or a single number) giving the
#'   significance level for bootstrap intervals. Defaults to 0.05.
#'
#' @return \item{values}{A vector of estimated eigenvalues \eqn{\hat \nu_1 ,
#'   \ldots , \hat \nu_M}.} \item{functions}{A
#'   \code{\link[funData]{multiFunData}} object containing the estimated
#'   multivariate functional principal components \eqn{\hat \psi_1, \ldots, \hat
#'   \psi_M}.} \item{scores}{ A matrix of dimension \code{N x M} containing the
#'   estimated scores \eqn{\hat \rho_{im}}.} \item{Yhat}{A
#'   \code{\link[funData]{multiFunData}} object containing estimated
#'   trajectories for each observation based on the truncated Karhunen-Lo\`{e}ve
#'   representation and the estimated scores and eigenfunctions.} \item{CI}{A
#'   list of the same length as \code{bootstrapAlpha}, containing the pointwise
#'   lower and upper bootstrap confidence bands for each significance level in
#'   form of \code{\link[funData]{multiFunData}} objects (only if
#'   \code{bootstrap = TRUE}).}
#'
#' @export MFPCA
#'
#' @importFrom foreach %do%
#'
#' @seealso \code{\link[funData]{multiFunData}}, \code{\link{PACE}},
#'   \code{\link{univBasisExpansion}}.
#'
#' @examples
#' oldPar <- par(no.readonly = TRUE)
#'
#' set.seed(1)
#'
#' ### simulate data (one-dimensional domains)
#' sim <-  simMultiFunData(type = "split", xVal = list(seq(0,1,0.01), seq(-0.5,0.5,0.02)),
#'                         M = 5, eFunType = "Poly", eValType = "linear", N = 100)
#'
#' # MFPCA based on univariate FPCA
#' \donttest{uFPCA <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"),
#'                                                                   list(type = "uFPCA")))}
#'
#' # MFPCA based on univariate spline expansions
#' splines <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "splines", k = 10),
#'                                                           list(type = "splines", k = 10)))
#'
#' # flip to make results more clear
#' \donttest{uFPCA$functions <- flipFuns(sim$trueFuns, uFPCA$functions)}
#' splines$functions <- flipFuns(sim$trueFuns, splines$functions)
#'
#' par(mfrow = c(1,2))
#' plot(sim$trueFuns[[1]], main = "Eigenfunctions", lwd = 2)
#' \donttest{plot(uFPCA$functions[[1]], lty = 2, add = TRUE)}
#' plot(splines$functions[[1]], lty = 3, add = TRUE)
#'
#' plot(sim$trueFuns[[2]], main = "Eigenfunctions", lwd = 2)
#' \donttest{plot(uFPCA$functions[[2]], lty = 2, add = TRUE)}
#' plot(splines$functions[[2]], lty = 3, add = TRUE)
#' legend("bottomleft", c("True", "uFPCA", "splines"), lty = 1:3, lwd = c(2,1,1))
#'
#' ### simulate data (two- and one-dimensional domains)
#' \donttest{
#' ### ATTENTION: Takes long
#' sim <-  simMultiFunData(type = "weighted",
#'                  xVal = list(list(seq(0,1,0.01), seq(-1,1,0.02)), list(seq(-0.5,0.5,0.01))),
#'                  M = list(c(4,5), 20), eFunType = list(c("Fourier", "Fourier"), "Poly"),
#'                  eValType = "exponential", N = 150)
#'
#' # MFPCA based on univariate spline expansions (for images) and univariate FPCA (for functions)
#' pca <- MFPCA(sim$simData, M = 10, uniExpansions = list(list(type = "splines", k = c(10,12)),
#'                                                        list(type = "uFPCA")))
#'
#' # flip to make results more clear
#' pca$functions <- flipFuns(extractObs(sim$trueFuns, obs = 1:10), pca$functions)
#'
#' par(mfrow = c(5,2), mar = rep(2,4))
#' for(m in 2:6) # for m = 1, image.plot (used in plot(funData)) produces an error...
#' {
#'   plot(sim$trueFuns[[1]], main = paste("True, m = ", m), obs = m)
#'   plot(pca$functions[[1]], main = paste("Estimate, m = ", m), obs = m)
#' }
#'
#' par(mfrow = c(1,1))
#' plot(sim$trueFuns[[2]], main = "Eigenfunctions", lwd = 2, obs=  1:5)
#' plot(pca$functions[[2]], lty = 2, add = TRUE, obs=  1:5)
#' legend("bottomleft", c("True", "MFPCA"), lty = 1:2, lwd = c(2,1))
#' }
#' par(oldPar)
MFPCA <- function(mFData, M, uniExpansions, weights = rep(1, length(mFData)), Yhat = FALSE, bootstrap = FALSE, nBootstrap = NULL, bootstrapAlpha = 0.05)
{
  # number of components
  p <- length(mFData)

  # number of observations
  N <- nObs(mFData)

  if(length(uniExpansions) != p)
    stop("Function MFPCA_multidim: multivariate functional data object and univariate expansions must have the same length!")

  if(bootstrap)
  {
    if(is.null(nBootstrap))
      stop("Specify number of bootstrap iterations")

    if(any(!(0 < bootstrapAlpha & bootstrapAlpha < 1)))
      stop("Significance level for bootstrap confidence bands must be in (0,1).")

  }

  # dimension for each component
  dimSupp <- dimSupp(mFData)

  # calculate univariate basis expansion for all components (if necessary)
  #   uniBasis <- vector("list", p)

  # for each component: find univariate basis expansion
  #   for(j in 1:p)
  #   {
  #     if(all(c("scores", "functions") %in% names(uniExpansions[[j]])))
  #     {
  #       tmp <- uniExpansions[[j]]
  #       type <- "user"
  #     }
  #     else
  #     {
  #       tmp <- switch(uniExpansions[[j]]$type,
  #                     "uFPCA" = do.call(PACE, args = findUniArgs(uniExpansions[[j]], mFData[[j]])),
  #                     "splines" = do.call(univBasisExpansion,  args = c(findUniArgs(uniExpansions[[j]], mFData[[j]]), pen = FALSE)),
  #                     "splinesPen" = do.call(univBasisExpansion,  args = c(findUniArgs(uniExpansions[[j]], mFData[[j]]), pen = TRUE)),
  #                     stop("Function MFPCA: uniExpansions type must be either 'uFPCA', 'splines' or 'splinesPen'")
  #       )
  #
  #       type <- uniExpansions[[j]]$type
  #     }
  #
  #     uniBasis[[j]] <- list(type = type, scores = tmp$scores, functions = tmp$functions@X)
  #
  #     if(dimSupp[j] == 2)
  #       uniBasis[[j]]$basisLong <- tmp$basisLong
  #   }

  uniBasis <- lapply(uniExpansion, function(l){univDecomp(type = l$type, data = l$data, params = l$params)})

  type <- sapply(uniExpansion, function(l){l$type})

  # Multivariate FPCA
  npc <- sapply(uniBasis, function(x){dim(x$scores)[2]}) # get number of univariate basis functions

  if(M > sum(npc))
    stop("Function MFPCA: total number of univariate basis functions must be greater or equal M!")

  #
  #   #  Block matrix of scalar products for each basis
  #   B <- array(0, dim = c(sum(npc), sum(npc)))
  #
  #   for(j in 1:p) # calculate block-wise
  #   {
  #     if(uniExpansions[[j]]$type == "uFPCA") # ONB -> matrix of scalar products is just the identity
  #       B[tmp[j]+ 1: npc[j], tmp[j] + 1:npc[j]] <- diag(npc[j])
  #     else # calculate scalar products
  #       B[tmp[j]+ 1: npc[j], tmp[j] + 1:npc[j]] <- .calcBasisIntegrals(uniBasis[[j]]$functions, dimSupp[j], mFData[[j]]@xVal)
  #   }

#   Z <- allScores %*% diag(allWeights)
#
#   # check if non-orthonormal basis functions used
   if(all(foreach::foreach(j = 1:p, .combine = "c")%do%{uniBasis$ortho}))
   {
     Bchol = NULL
#     tmpSVD <- irlba::irlba(1/sqrt(N-1) * Z, nv = M)
#
#     vectors <- tmpSVD$v
#     values <- tmpSVD$d
   }
   else
   {
#     # Cholesky decomposition of B = block diagonal of Cholesky decompositions
     Bchol <- bdiag(lapply(uniBasis), function(l){ifelse(l$ortho, Diagonal(n = ncol(l$scores)), chol(l$B)})
#
#     tmpSVD <- irlba::irlba(1/sqrt(N-1) * Z %*% t(Bchol), nv = M)
#
#     vectors <- t(Bchol) %*% tmpSVD$v
#     values <- tmpSVD$d
   }
#
#   # normalization factors
#   normFactors <- 1/sqrt(diag(t(vectors) %*% crossprod(Z) %*% vectors)/(N-1))
#
#   # calculate scores
#   scores <- Z %*% vectors
#   scores <- scores %*% diag(sqrt(values) * normFactors) # normalization
#
#   # calculate eigenfunctions (incl. normalization)
#   tmpWeights <- 1/(N-1) *  crossprod(Z) %*% vectors
#   eFunctions <- foreach::foreach(j = 1:p){
#     univExpansion(type = uniExpansion[[j]]$type,
#                   scores = weights[j] * tmpWeights[npcCum[j]+1:npc[j],] %*%  diag(1/sqrt(values) * normFactors),
#                   xVal = mFData[[j]]$xVal,
#                   functions = uniBasis[[j]]$functions,
#                   params = uniBasis[[j]]$settings)
#   }
#
#   # calculate truncated Karhunen-Loeve representation
#   Yhat <- foreach::foreach(j = 1:p){
#     univExpansion(type = uniExpansion[[j]]$type,
#                   scores = scores[npcCum[j]+1:npc[j],],
#                   xVal = mFData[[j]]$xVal,
#                   functions = uniBasis[[j]]$functions,
#                   params = uniBasis[[j]]$settings)
#   }


  res <- calcMFPCA(N = N, p = p, Bchol = Bchol, M = M, type = type, weights = weights,
                   npc = npc, xVal = getxVal(mFData), uniBasis = uniVasis, Yhat = Yhat)

  #   # and then calculate covariance for each combination, component-wise multiplication with weights
  #   Z <- cov(allScores) * sqrt(outer(allWeights, allWeights, "*"))
  #
  #   # do eigendecomposition
  #   C <- eigen(B%*%Z)

  # factors for normalizations
#   normFactors <- diag(t(C$vectors) %*% Z %*% C$vectors)[1:M]
#
#   # calculate (multivariate) scores
#   scores <- Re(allScores %*% diag(sqrt(allWeights)) %*% C$vectors[, 1:M] %*% diag(sqrt(C$values[1:M])/sqrt(normFactors)) )
#
#   # calculate multivariate eigenfunctions and truncated Karhunen-Lo\`{e}ve representation
#   eFunctions <- vector("list", p)
#   Yhat <- vector("list", p)
#
#   for(j in 1:p)
#   {
#     # calculate eigenfunctions
#     if(dimSupp[j] == 1) # one-dimensional function
#       eFuns <-  1/sqrt(weights[j]) * t(Z[tmp[j]+1:npc[j], ] %*% C$vectors[, 1:M]) %*% uniBasis[[j]]$functions
#     else # two-dimesional function (otherwise function stops before!)
#       eFuns <-  1/sqrt(weights[j]) * t(uniBasis[[j]]$basisLong %*% Z[tmp[j]+1:npc[j], ] %*% C$vectors[, 1:M])
#
#     # normalize
#     eFuns <- Re(diag(1/sqrt(C$values[1:M] * normFactors)) %*% eFuns)
#
#     # truncated Karhunen-Loève representation (reconstruction)
#     recons <-  scores %*%  eFuns
#
#     # for two-dimensional functions: reshape eigenfunctions and reconstruction to image
#     if(dimSupp[j] == 2) # two-dimensional function: reshape
#     {
#       eFuns <- array(eFuns, dim = c(M,length(mFData[[j]]@xVal[[1]]), length(mFData[[j]]@xVal[[2]])))
#       recons <-  array(recons, dim = c(N, length(mFData[[j]]@xVal[[1]]), length(mFData[[j]]@xVal[[2]])))
#     }
#
#     eFunctions[[j]] <- funData(xVal = mFData[[j]]@xVal, X = eFuns)
#     Yhat[[j]] <- funData(xVal = mFData[[j]]@xVal, X = recons)
#   }




  # bootstrap for eigenfunctions
  if(bootstrap)
  {
    booteFuns <- vector("list", p)

    for(j in 1:p)
      booteFuns[[j]] <- array(NA, dim  = c(nBootstrap, M, sapply(mFData[[j]]@xVal, length)))

    for(n in 1:nBootstrap)
    {
      bootObs <- sample(N, replace = TRUE)

      bootBasis <- vector("list", p)

      for(j in 1:p)
      {
        if(uniExpansion[[j]]$type == "uFPCA") # re-estimate scores AND functions
          bootBasis[[j]] <- univDecomp(type = uniExpansion[[j]]$type, data = extractObs(mFData, obs = bootObs), params = uniExpansion[[j]]$params)
        else # resample scores (functions are given and scores can simply be resampled)
          bootBasis[[j]] <- list(scores = uniBasis[[j]]$scores[bootObs, ], B = uniBasis[[j]]$B, ortho = uniBasis[[j]]$ortho, functions = uniBasis[[j]]$functions)
      }

      npcBoot <- sapply(bootBasis, function(x){dim(x$scores)[2]}) # get number of univariate basis functions

      if(M > sum(npcBoot))
        stop("Function MFPCA (bootstrap): total number of univariate basis functions must be greater or equal M!")

#       ### calculate MFPCA -> the same as above for uniBasis now for bootBasis (no Yhat)
#
#       # combine all scores of bootstrap observations
#       bootScores <- foreach::foreach(j = 1:p, .combine = "cbind")%do%{bootBasis[[j]]$scores}
#
#       Z <- cov(bootScores) * sqrt(outer(allWeights, allWeights, "*")) # and then calculate covariance for each combination, component-wise multiplication with weights
#
#       # do eigendecomposition
#       C <- eigen(B%*%Z)
#
#       # factors for normalizations
#       normFactors <- diag(t(C$vectors) %*% Z %*% C$vectors)[1:M]
#
#       # calculate multivariate eigenfunctions for bootstrap sample
#
#       for(j in 1:p)
#       {
#         # calculate eigenfunctions
#         if(dimSupp[j] == 1) # one-dimensional function
#           eFuns <-  1/sqrt(weights[j]) * t(Z[tmp[j]+1:npc[j], ] %*% C$vectors[, 1:M]) %*% bootBasis[[j]]$functions
#         else # two-dimensional function (otherwise function stops before!)
#           eFuns <-  1/sqrt(weights[j]) * t(bootBasis[[j]]$basisLong %*% Z[tmp[j]+1:npc[j], ] %*% C$vectors[, 1:M])
#
#         # normalize
#         eFuns <- Re(diag(1/sqrt(C$values[1:M] * normFactors)) %*% eFuns)
#
#         # save in temporary list (must check for flipping before final save)
#         if(dimSupp[j] == 1)
#           tmpFuns[[j]] <- funData(mFData[[j]]@xVal, eFuns)
#         else # two-dimensional function: reshape
#           tmpFuns[[j]] <- funData(mFData[[j]]@xVal,
#                                   array(eFuns, dim = c(M,length(mFData[[j]]@xVal[[1]]), length(mFData[[j]]@xVal[[2]]))))
#
#       }

    # calculate MFPCA for bootstrap sample (Bchol must not be recalculated as uFPCA basis functions are orthonormal!)
    tmpFuns <- calcMFPCA(N = N, p = p, Bchol = Bchol, M = M, type = type, weights = weights,
                         npc = npcBoot, xVal = getxVal(mFData), uniBasis = bootBasis, Yhat = FALSE)$functions

      # flip bootstrap estimates if necessary
      tmpFuns <- flipFuns(res$functions, tmpFuns)

      # save in booteFuns
      for(j in 1:p)
      {
        if(dimSupp[j] == 1)
          booteFuns[[j]][n,,] <- tmpFuns[[j]]@X
        else # two-dimensional function
          booteFuns[[j]][n,,,] <- tmpFuns[[j]]@X
      }
    }

    CI <- vector("list", length(bootstrapAlpha))

    for(alpha in 1:length(bootstrapAlpha))
    {
      bootCI_lower <- bootCI_upper <-  vector("list", p)
      for(j in 1:p)
      {
        bootCI_lower[[j]] <- funData(mFData[[j]]@xVal, apply(booteFuns[[j]], 2:length(dim(booteFuns[[j]])),
                                                             quantile, bootstrapAlpha[alpha]/2))
        bootCI_upper[[j]] <- funData(mFData[[j]]@xVal, apply(booteFuns[[j]],  2:length(dim(booteFuns[[j]])),
                                                             quantile, 1 - bootstrapAlpha[alpha]/2))
      }

      CI[[alpha]]$lower <- multiFunData(bootCI_lower)
      CI[[alpha]]$upper <- multiFunData(bootCI_upper)

      names(CI)[alpha] <- paste("alpha", bootstrapAlpha[alpha], sep = "_")
    }

    res$CI <- CI
  }

  return(res)
}

#' Internal function that implements the MFPCA algorithm for given univariate decompositions
#'
#' @keywords internal
calcMFPCA <- function(N, p, Bchol, M, type, weights, npc, xVal, uniBasis, Yhat = FALSE)
{
  # combine all scores
  allScores <- foreach::foreach(j = 1:p, .combine = "cBind")%do%{uniBasis[[j]]$scores}

  # de-mean scores (column-wise)
  allScores <- apply(allScores, 2, function(x){x - mean(x)})

  # block vector of weights
  allWeights <- foreach::foreach(j = 1:p, .combine = "c")%do%{rep(weights[j], npc[j])}

  Z <- allScores %*% diag(allWeights)

  # check if non-orthonormal basis functions used and calculate PCA on scores
  if(is.null(Bchol))
  {
    tmpSVD <- irlba::irlba(1/sqrt(N-1) * Z, nv = M)

    vectors <- tmpSVD$v
    values <- tmpSVD$d
  }
  else
  {
    tmpSVD <- irlba::irlba(1/sqrt(N-1) * Z %*% t(Bchol), nv = M)

    vectors <- t(Bchol) %*% tmpSVD$v
    values <- tmpSVD$d
  }

  # normalization factors
  normFactors <- 1/sqrt(diag(t(vectors) %*% crossprod(Z) %*% vectors)/(N-1))

  ### Calculate scores
  scores <- Z %*% vectors
  scores <- scores %*% diag(sqrt(values) * normFactors) # normalization

  ### Calculate eigenfunctions (incl. normalization)
  npcCum <- cumsum(c(0, npc)) # indices for blocks (-1)

  tmpWeights <- 1/(N-1) *  crossprod(Z) %*% vectors
  eFunctions <- foreach::foreach(j = 1:p){
    univExpansion(type = type[j],
                  scores = weights[j] * tmpWeights[npcCum[j]+1:npc[j],] %*%  diag(1/sqrt(values) * normFactors),
                  xVal = xVal[[j]],
                  functions = uniBasis[[j]]$functions,
                  params = uniBasis[[j]]$settings)
  }

  res <- list(values = values,
              functions = multiFunData(eFunctions),
              scores = scores)

  if(Yhat)
  {
    # calculate truncated Karhunen-Loeve representation
    Yhat <- foreach::foreach(j = 1:p){
      univExpansion(type = type[j],
                    scores = scores[npcCum[j]+1:npc[j],],
                    xVal = xVal[[j]],
                    functions = uniBasis[[j]]$functions,
                    params = uniBasis[[j]]$settings)
    }
    res$Yhat <- multiFunData(Yhat)
  }

  return(res)
}


#' Utility function for defining arguments of univariate basis expansion
#'
#' This is a utility function that defines the arguments
#' of\code{\link{univBasisExpansion}} in \code{\link{MFPCA}}. This is used as an
#' internal function in \code{\link{MFPCA}}.
#'
#' @param uniExpansion A list corresponding to one univariate expansion as
#'   described in \code{\link{MFPCA}}.
#' @param funDataObject An object of class \code{\link[funData]{funData}},
#'   representing the (univariate) functional data, for which the expansion will
#'   be calculated
#'
#' @return \item{args}{A list of parameters to be passed to
#'   \code{\link{univBasisExpansion}}}.
#'
#' @keywords internal
findUniArgs <- function(uniExpansion, funDataObject)
{
  if(all(names(uniExpansion) == "type")) # i.e. all other values are defaults
    args <- list(funDataObject = funDataObject, NULL)
  else
  {
    args <- uniExpansion[names(uniExpansion) != "type"]
    args$funDataObject <- funDataObject
  }

  return(args)
}
