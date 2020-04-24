#' @import funData
NULL

# define global variable j, used by the foreach package and confusing R CMD CHECK
globalVariables('j')

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
#' @param argvals List of corresponding x-values.
#'
#' @return A matrix containing the scalar product of all combinations of basis functions (matrix \eqn{B^{(j)}})
#'
#' @seealso \code{\link{MFPCA}}, \code{\link[funData]{dimSupp}}
#'
#' @keywords internal
calcBasisIntegrals <- function(basisFunctions, dimSupp, argvals)
{
  npc <- dim(basisFunctions)[1]

  #  integral basis matrix
  B <- array(0, dim = c(npc, npc))


  if(dimSupp == 1) # one-dimensional domain
  {
    w <- funData::.intWeights(argvals[[1]])

    for(m in seq_len(npc))
    {
      for(n in seq_len(m))
        B[m, n] <- B[n, m] <- (basisFunctions[m, ]* basisFunctions[n, ])%*%w
    }
  }
  else # two-dimesional domain (otherwise function stops before!)
  {
    w1 <- t(funData::.intWeights(argvals[[1]]))
    w2 <- funData::.intWeights(argvals[[2]])

    for(m in seq_len(npc))
    {
      for(n in seq_len(m))
        B[m, n] <- B[n, m] <-  w1 %*%(basisFunctions[m, , ]* basisFunctions[n, ,])%*%w2
    }
  }

  return(B)
}


#' Multivariate functional principal component analysis for functions on
#' different (dimensional) domains
#'
#' This function calculates a multivariate functional principal component
#' analysis (MFPCA) based on i.i.d. observations \eqn{x_1, \ldots, x_N} of
#' a multivariate functional data-generating process \eqn{X = (X^{(1)},
#' \ldots X^{(p)})}{X = X^(1), \ldots, X^(p)} with elements \eqn{X^{(j)}
#' \in L^2(\mathcal{T}_j)}{X^(j) in L^2(calT_j)} defined on a domain
#' \eqn{\mathcal{T}_j \subset IR^{d_j}}{calT_j of IR^{d_j}}. In
#' particular, the elements can be defined on different (dimensional)
#' domains. The results contain the mean function, the estimated
#' multivariate functional principal components \eqn{\hat \psi_1, \ldots,
#' \hat \psi_M} (having the same structure as \eqn{x_i}), the associated
#' eigenvalues \eqn{\hat \nu_1 \geq \ldots \geq \hat \nu_M > 0} and the
#' individual scores \eqn{\hat \rho_{im} = \widehat{<x_i, \psi_m>}}{\hat
#' \rho_{im} = \hat{<x_i, \psi_m>}}. Moreover, estimated trajectories for
#' each observation based on the truncated Karhunen-Loeve representation
#' \deqn{\hat x_i = \sum_{m = 1}^M \hat \rho_{im} \hat \psi_m}{\hat x_i =
#' \sum_{m = 1}^M \hat \rho_{im} \hat \psi_m} are given if desired
#' (\code{fit = TRUE}). The implementation of the observations \eqn{x_i =
#' (x_i^{(1)}, \ldots , x_i^{(p)}),~ i = 1 , \ldots, N}{x_i = (x_i^(1),
#' \ldots , x_i^(p)), i = 1 , \ldots, N}, the mean function and
#' multivariate functional principal components \eqn{\hat \psi_1, \ldots,
#' \hat \psi_M} uses the \code{\link[funData]{multiFunData}} class, which
#' is defined in the package \pkg{funData}.
#'
#' \subsection{Weighted MFPCA}{If the elements vary considerably in
#' domain, range or variation, a weight vector \eqn{w_1 , \ldots, w_p} can
#' be supplied and the MFPCA is based on the weighted scalar product
#' \deqn{<<f,g>>_w = \sum_{j = 1}^p w_j \int_{\mathcal{T}_j} f^{(j)}(t)
#' g^{(j)}(t) \mathrm{d} t}{<<f,g>>_w = \sum_{j = 1}^p w_j \int_\calT_j
#' f^(j)(t) g^(j)(t) d t} and the corresponding weighted covariance
#' operator \eqn{\Gamma_w}.}
#'
#' \subsection{Bootstrap}{If \code{bootstrap = TRUE}, pointwise bootstrap
#' confidence bands are generated for the multivariate eigenvalues
#' \eqn{\hat \nu_1, \ldots, \hat \nu_M } as well as for multivariate
#' functional principal components \eqn{\hat \psi_1, \ldots, \hat
#' \psi_M}{\hat \psi_1, \ldots, \hat \psi_M}. The parameter
#' \code{nBootstrap} gives the number of bootstrap iterations. In each
#' iteration, the observations are resampled on the level of
#' (multivariate) functions and the whole MFPCA is recalculated. In
#' particular, if the univariate basis depends on the data (FPCA
#' approaches), basis functions and scores are both re-estimated. If the
#' basis functions are fixed (e.g. splines), the scores from the original
#' estimate are used to speed up the calculations. The confidence bands
#' for the eigenfunctions are calculated separately for each element as
#' pointwise percentile bootstrap confidence intervals. Analogously, the
#' confidence bands for the eigenvalues are also percentile bootstrap
#' confidence bands. The significance level(s) can be defined by the
#' \code{bootstrapAlpha} parameter, which defaults to 5\%. As a result,
#' the \code{MFPCA} function returns a list \code{CI} of the same length
#' as \code{bootstrapAlpha}, containing the lower and upper bounds of the
#' confidence bands for the principal components as \code{multiFunData}
#' objects of the same structure as \code{mFData}. The confidence bands
#' for the eigenvalues are returned in a list \code{CIvalues}, containing
#' the upper and lower bounds for each significance level.}
#'
#'
#' \subsection{Univariate Expansions}{The multivariate functional
#' principal component analysis relies on a univariate basis expansion for
#' each element \eqn{X^{(j)}}{X^(j)}. The univariate basis representation
#' is calculated using the \code{\link{univDecomp}} function, that passes
#' the univariate functional observations and optional parameters to the
#' specific function. The univariate decompositions are specified via the
#' \code{uniExpansions} argument in the \code{MFPCA} function. It is a
#' list of the same length as the \code{mFData} object, i.e. having one
#' entry for each element of the multivariate functional data. For each
#' element, \code{uniExpansion} must specify at least the type of basis
#' functions to use. Additionally, one may add further parameters. The
#' following basis representations are supported: \itemize{\item Given
#' basis functions. Then \code{uniExpansions[[j]] = list(type = "given",
#' functions, scores, ortho)}, where \code{functions} is a \code{funData}
#' object on the same domain as \code{mFData}, containing the given basis
#' functions. The parameters \code{scores} and \code{ortho} are optional.
#' \code{scores} is an \code{N x K} matrix containing the scores (or
#' coefficients) of the observed functions for the given basis functions,
#' where \code{N} is the number of observed functions and \code{K} is the
#' number of basis functions. If this is not supplied, the scores are
#' calculated. The parameter \code{ortho} specifies whether the given
#' basis functions are orthonormal \code{orhto = TRUE} or not \code{ortho
#' = FALSE}. If \code{ortho} is not supplied, the functions are treated as
#' non-orthogonal. \code{scores} and \code{ortho} are not checked for
#' plausibility, use them on your own risk!  \item Univariate functional
#' principal component analysis. Then \code{uniExpansions[[j]] = list(type
#' = "uFPCA", nbasis, pve, npc, makePD)}, where
#' \code{nbasis,pve,npc,makePD} are parameters passed to the
#' \code{\link{PACE}} function for calculating the univariate functional
#' principal component analysis. \item Basis functions expansions from the
#' package \pkg{fda}. Then \code{uniExpansions[[j]] = list(type = "fda",
#' ...)}, where \code{...} are passed to
#' \code{\link[funData]{funData2fd}}, which heavily builds on
#' \code{\link[fda]{eval.fd}}. If \pkg{fda} is not available, a warning is
#' thrown. \item Spline basis functions (not penalized). Then
#' \code{uniExpansions[[j]] = list(type = "splines1D", bs, m, k)}, where
#' \code{bs,m,k} are passed to the functions \code{\link{univDecomp}} and
#' \code{\link{univExpansion}}. For two-dimensional tensor product
#' splines, use \code{type = "splines2D"}. \item Spline basis functions
#' (with smoothness penalty). Then \code{uniExpansions[[j]] = list(type =
#' "splines1Dpen", bs, m, k)}, where \code{bs,m,k} are passed to the
#' functions \code{\link{univDecomp}} and \code{\link{univExpansion}}.
#' Analogously to the unpenalized case, use \code{type = "splines2Dpen"}
#' for 2D penalized tensor product splines. \item Cosine basis functions.
#' Use \code{uniExpansions[[j]] = list(type = "DCT2D", qThresh, parallel)}
#' for functions one two-dimensional domains (images) and \code{type =
#' "DCT3D"} for 3D images. The calculation is based on the discrete cosine
#' transform (DCT) implemented in the C-library \code{fftw3}. If this
#' library is not available, the function will throw  a warning.
#' \code{qThresh} gives the quantile for hard thresholding the basis
#' coefficients based on their absolute value. If \code{parallel = TRUE},
#' the coefficients for different images are calculated in parallel.} See
#' \code{\link{univDecomp}} and \code{\link{univExpansion}} for details.}
#'
#' @param mFData A  \code{\link[funData]{multiFunData}} object containing
#'   the \code{N} observations.
#' @param M The number of multivariate functional principal components to
#'   calculate.
#' @param uniExpansions A list characterizing the (univariate) expansion
#'   that is calculated for each element. See Details.
#' @param weights An optional vector of weights, defaults to \code{1} for
#'   each element. See Details.
#' @param fit Logical. If \code{TRUE}, a truncated multivariate
#'   Karhunen-Loeve representation for the data is calculated based on the
#'   estimated scores and eigenfunctions.
#' @param approx.eigen Logical. If \code{TRUE}, the eigenanalysis problem
#'   for the estimated covariance matrix is solved approximately using the
#'   \pkg{irlba} package, which is much faster. If the number \code{M} of
#'   eigenvalues to calculate is high with respect to the number of
#'   observations in \code{mFData} or the number of estimated univariate
#'   eigenfunctions, the approximation may be inappropriate. In this case,
#'   approx.eigen is set to \code{FALSE} and the function throws a
#'   warning. Defaults to \code{FALSE}.
#' @param bootstrap Logical. If \code{TRUE}, pointwise bootstrap
#'   confidence bands are calculated for the multivariate functional
#'   principal components. Defaults to \code{FALSE}. See Details.
#' @param nBootstrap The number of bootstrap iterations to use. Defaults
#'   to \code{NULL}, which leads to an error, if \code{bootstrap = TRUE}.
#' @param bootstrapAlpha A vector of numerics (or a single number) giving
#'   the significance level for bootstrap intervals. Defaults to
#'   \code{0.05}.
#' @param bootstrapStrat A stratification variable for bootstrap. Must be
#'   a factor of length \code{nObs(mFData)} or \code{NULL} (default). If
#'   \code{NULL}, no stratification is made in the bootstrap resampling,
#'   i.e. the curves are sampled with replacement. If
#'   \code{bootstrapStrat} is not \code{NULL}, the curves are resampled
#'   with replacement within the groups defined by \code{bootstrapStrat},
#'   hence keeping the group proportions fixed.
#' @param verbose Logical. If \code{TRUE}, the function reports
#'   extra-information about the progress (incl. timestamps). Defaults to
#'   \code{options()$verbose}.
#'
#' @return An object of class \code{MFPCAfit} containing the following
#'   components: \item{values}{A vector of estimated eigenvalues \eqn{\hat
#'   \nu_1 , \ldots , \hat \nu_M}.} \item{functions}{A
#'   \code{\link[funData]{multiFunData}} object containing the estimated
#'   multivariate functional principal components \eqn{\hat \psi_1,
#'   \ldots, \hat \psi_M}.} \item{scores}{ A matrix of dimension \code{N x
#'   M} containing the estimated scores \eqn{\hat \rho_{im}}.}
#'   \item{vectors}{A matrix representing the eigenvectors associated with
#'   the combined univariate score vectors. This might be helpful for
#'   calculating predictions.} \item{normFactors}{The normalizing factors
#'   used for calculating the multivariate eigenfunctions and scores. This
#'   might be helpful when calculation predictions.} \item{meanFunction}{A
#'   multivariate functional data object, corresponding to the mean
#'   function. The MFPCA is applied to the de-meaned functions in
#'   \code{mFData}.}\item{fit}{A \code{\link[funData]{multiFunData}}
#'   object containing estimated trajectories for each observation based
#'   on the truncated Karhunen-Loeve representation and the estimated
#'   scores and eigenfunctions.} \item{CI}{A list of the same length as
#'   \code{bootstrapAlpha}, containing the pointwise lower and upper
#'   bootstrap confidence bands for each eigenfunction and each
#'   significance level in form of \code{\link[funData]{multiFunData}}
#'   objects (only if \code{bootstrap = TRUE}).} \item{CIvalues}{A list of
#'   the same length as \code{bootstrapAlpha}, containing the lower and
#'   upper bootstrap confidence bands for each eigenvalue and each
#'   significance level (only if \code{bootstrap = TRUE}).}
#'
#' @export MFPCA
#'
#' @importFrom foreach %do%
#' @importFrom utils packageVersion
#'
#' @references C. Happ, S. Greven (2018): Multivariate Functional
#'   Principal Component Analysis for Data Observed on Different
#'   (Dimensional) Domains. Journal of the American Statistical
#'   Association, 113(522): 649-659. DOI:
#'   \doi{10.1080/01621459.2016.1273115}
#'
#' @references C. Happ-Kurz (2020): Object-Oriented Software for Functional
#'   Data. Journal of Statistical Software, 93(5): 1-38. DOI:
#'   \doi{10.18637/jss.v093.i05}
#'
#' @seealso See Happ-Kurz (2020. \doi{10.18637/jss.v093.i05}) for a general
#'   introduction to the \pkg{funData} package and it's interplay with
#'   \pkg{MFPCA}. This file also includes a case study on how to use
#'   \code{MFPCA}. Useful functions: \code{\link[funData]{multiFunData}},
#'   \code{\link{PACE}}, \code{\link{univDecomp}},
#'   \code{\link{univExpansion}},
#'   \code{\link[=summary.MFPCAfit]{summary}},
#'   \code{\link[=plot.MFPCAfit]{plot}},
#'   \code{\link[=scoreplot.MFPCAfit]{scoreplot}}
#'
#' @examples
#' oldPar <- par(no.readonly = TRUE)
#'
#' set.seed(1)
#'
#' ### simulate data (one-dimensional domains)
#' sim <-  simMultiFunData(type = "split", argvals = list(seq(0,1,0.01), seq(-0.5,0.5,0.02)),
#'                         M = 5, eFunType = "Poly", eValType = "linear", N = 100)
#'
#' # MFPCA based on univariate FPCA
#' uFPCA <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"),
#'                                                                   list(type = "uFPCA")))

#' summary(uFPCA)
#' plot(uFPCA) # plot the eigenfunctions as perturbations of the mean
#' scoreplot(uFPCA) # plot the scores
#' 
#' # MFPCA based on univariate spline expansions
#' splines <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "splines1D", k = 10),
#'                                                           list(type = "splines1D", k = 10)),
#'                  fit = TRUE) # calculate reconstruction, too
#' summary(splines)
#' plot(splines) # plot the eigenfunctions as perturbations of the mean
#' scoreplot(splines) # plot the scores
#' 
#' ### Compare estimates to true eigenfunctions
#' # flip to make results more clear
#' uFPCA$functions <- flipFuns(sim$trueFuns, uFPCA$functions)
#' splines$functions <- flipFuns(sim$trueFuns, splines$functions)
#' 
#' par(mfrow = c(1,2))
#' plot(sim$trueFuns[[1]], main = "Eigenfunctions\n1st Element", lwd = 2)
#' plot(uFPCA$functions[[1]], lty = 2, add = TRUE)
#' plot(splines$functions[[1]], lty = 3, add = TRUE)
#' 
#' plot(sim$trueFuns[[2]], main = "Eigenfunctions\n2nd Element", lwd = 2)
#' plot(uFPCA$functions[[2]], lty = 2, add = TRUE)
#' plot(splines$functions[[2]], lty = 3, add = TRUE)
#' legend("bottomleft", c("True", "uFPCA", "splines"), lty = 1:3, lwd = c(2,1,1))
#' 
#' # Test reconstruction for the first 10 observations
#' plot(sim$simData[[1]], obs = 1:10, main = "Reconstruction\n1st Element", lwd = 2)
#' plot(splines$fit[[1]], obs = 1:10, lty = 2, col = 1, add = TRUE)
#' 
#' plot(sim$simData[[2]], obs = 1:10, main = "Reconstruction\n2nd Element", lwd = 2)
#' plot(splines$fit[[2]], obs = 1:10, lty = 2, col = 1, add = TRUE)
#' legend("bottomleft", c("True", "Reconstruction"), lty = c(1,2), lwd = c(2,1))
#' 
#' # MFPCA with Bootstrap-CI for the first 2 eigenfunctions
#' ### ATTENTION: Takes long
#' \donttest{
#' splinesBoot <- MFPCA(sim$simData, M = 2, uniExpansions = list(list(type = "splines1D", k = 10),
#'                                                           list(type = "splines1D", k = 10)),
#'                  bootstrap = TRUE, nBootstrap = 100, bootstrapAlpha = c(0.05, 0.1), verbose = TRUE)
#' summary(splinesBoot)
#'                                  
#' plot(splinesBoot$functions[[1]], ylim = c(-2,1.5))
#' plot(splinesBoot$CI$alpha_0.05$lower[[1]], lty = 2, add = TRUE)
#' plot(splinesBoot$CI$alpha_0.05$upper[[1]], lty = 2, add = TRUE)
#' plot(splinesBoot$CI$alpha_0.1$lower[[1]], lty = 3, add = TRUE)
#' plot(splinesBoot$CI$alpha_0.1$upper[[1]], lty = 3, add = TRUE)
#' abline(h = 0, col = "gray")
#'  
#' plot(splinesBoot$functions[[2]], ylim = c(-1,2.5))
#' plot(splinesBoot$CI$alpha_0.05$lower[[2]], lty = 2, add = TRUE)
#' plot(splinesBoot$CI$alpha_0.05$upper[[2]], lty = 2, add = TRUE)
#' plot(splinesBoot$CI$alpha_0.1$lower[[2]], lty = 3, add = TRUE)
#' plot(splinesBoot$CI$alpha_0.1$upper[[2]], lty = 3, add = TRUE)
#' abline(h = 0, col = "gray")
#' legend("topleft", c("Estimate", "95% CI", "90% CI"), lty = 1:3, lwd = c(2,1,1))
#' 
#' # Plot 95% confidence bands for eigenvalues
#' plot(1:2, splinesBoot$values, pch = 20, ylim = c(0, 1.5), 
#'      main = "Estimated eigenvalues with 95% CI",
#'      xlab = "Eigenvalue no.", ylab = "")
#' arrows(1:2, splinesBoot$CIvalues$alpha_0.05$lower,
#'        1:2, splinesBoot$CIvalues$alpha_0.05$upper,
#'        length = 0.05, angle = 90, code = 3)
#' points(1:2, sim$trueVals[1:2], pch = 20, col = 4)
#' legend("topright", c("Estimate", "True value"), pch = 20, col = c(1,4))
#' }
#' 
#' ### simulate data (two- and one-dimensional domains)
#' ### ATTENTION: Takes long
#' \donttest{
#' set.seed(2)
#' sim <-  simMultiFunData(type = "weighted",
#'                  argvals = list(list(seq(0,1,0.01), seq(-1,1,0.02)), list(seq(-0.5,0.5,0.01))),
#'                  M = list(c(4,5), 20), eFunType = list(c("Fourier", "Fourier"), "Poly"),
#'                  eValType = "exponential", N = 150)
#' 
#' # MFPCA based on univariate spline expansions (for images) and univariate FPCA (for functions)
#' pca <- MFPCA(sim$simData, M = 10,
#'              uniExpansions = list(list(type = "splines2D", k = c(10,12)),
#'                              list(type = "uFPCA")))
#' summary(pca)
#' plot(pca) # plot the eigenfunctions as perturbations of the mean
#' scoreplot(pca) # plot the scores
#' 
#' ### Compare to true eigenfunctions
#' # flip to make results more clear
#' pca$functions <- flipFuns(sim$trueFuns[1:10], pca$functions)
#' 
#' par(mfrow = c(5,2), mar = rep(2,4))
#' for(m in 2:6) # for m = 1, image.plot (used in plot(funData)) produces an error...
#' {
#'   plot(sim$trueFuns[[1]], main = paste("True, m = ", m), obs = m)
#'   plot(pca$functions[[1]], main = paste("Estimate, m = ", m), obs = m)
#' }
#' 
#' par(mfrow = c(1,1))
#' plot(sim$trueFuns[[2]], main = "Eigenfunctions (2nd element)", lwd = 2, obs=  1:5)
#' plot(pca$functions[[2]], lty = 2, add = TRUE, obs=  1:5)
#' legend("bottomleft", c("True", "MFPCA"), lty = 1:2, lwd = c(2,1))
#' }
#' par(oldPar)
MFPCA <- function(mFData, M, uniExpansions, weights = rep(1, length(mFData)), fit = FALSE, approx.eigen = FALSE,
                  bootstrap = FALSE, nBootstrap = NULL, bootstrapAlpha = 0.05, bootstrapStrat = NULL, 
                  verbose = options()$verbose)
{
  if(! inherits(mFData, "multiFunData"))
    stop("Parameter 'mFData' must be passed as a multiFunData object.")
  
  # number of components
  p <- length(mFData)
  # number of observations
  N <- nObs(mFData)
  
  if(!all(is.numeric(M), length(M) == 1, M > 0))
    stop("Parameter 'M' must be passed as a number > 0.")
  
  if(!(is.list(uniExpansions) & length(uniExpansions) == p))
    stop("Parameter 'uniExpansions' must be passed as a list with the same length as 'mFData'.")
  
  if(!(is.numeric(weights) & length(weights) == p))
    stop("Parameter 'weights' must be passed as a vector with the same length as 'mFData'.")
  
  if(!is.logical(fit))
    stop("Parameter 'fit' must be passed as a logical.")
  
  if(!is.logical(approx.eigen))
    stop("Parameter 'approx.eigen' must be passed as a logical.")
  
  if(!is.logical(bootstrap))
    stop("Parameter 'bootstrap' must be passed as a logical.")
  
  if(bootstrap)
  {
    if(is.null(nBootstrap))
      stop("Specify number of bootstrap iterations.")

    if(any(!(0 < bootstrapAlpha & bootstrapAlpha < 1)))
      stop("Significance level for bootstrap confidence bands must be in (0,1).")

    if(!is.null(bootstrapStrat))
    {
      if(!is.factor(bootstrapStrat))
        stop("bootstrapStrat must be either NULL or a factor.")
      
      if(length(bootstrapStrat) != nObs(mFData))
        stop("bootstrapStrat must have the same length as the number of observations in the mFData object.")
    }
  }
  
  if(!is.logical(verbose))
    stop("Parameter 'verbose' must be passed as a logical.")
  
  # dimension for each component
  dimSupp <- dimSupp(mFData)
  
  # get type of univariate expansions
  type <- vapply(uniExpansions, function(l){l$type}, FUN.VALUE = "")

  # de-mean functions -> coefficients are also de-meaned!
  # do not de-mean in uFPCA, as PACE gives a smooth estimate of the mean (see below)
  m <- meanFunction(mFData, na.rm = TRUE) # ignore NAs in data
  for(j in seq_len(p))
  { 
    if(type[j] != "uFPCA")
      mFData[[j]] <- mFData[[j]] - m[[j]]
  }

  if(verbose)
    cat("Calculating univariate basis expansions (", format(Sys.time(), "%T"), ")\n", sep = "")

  # calculate univariate basis expansion for all components
  uniBasis <- mapply(function(expansion, data){do.call(univDecomp, c(list(funDataObject = data), expansion))},
                     expansion = uniExpansions, data = mFData, SIMPLIFY = FALSE)

  # for uFPCA: replace estimated mean in m
  for(j in seq_len(p))
  {
    if(type[j] == "uFPCA")
      m[[j]] <- uniBasis[[j]]$meanFunction
  }

  # Multivariate FPCA
  npc <- vapply(uniBasis, function(x){dim(x$scores)[2]}, FUN.VALUE = 0) # get number of univariate basis functions

  if(M > sum(npc))
  {
    M <- sum(npc)
    warning("Function MFPCA: total number of univariate basis functions is smaller than given M. M was set to ", sum(npc), ".")
  } 

  # check if non-orthonormal basis functions used
  if(all(foreach::foreach(j = seq_len(p), .combine = "c")%do%{uniBasis[[j]]$ortho}))
    Bchol = NULL
  else
  {
    # Cholesky decomposition of B = block diagonal of Cholesky decompositions
    Bchol <- Matrix::bdiag(lapply(uniBasis, function(l){
      if(l$ortho)
        res <- Matrix::Diagonal(n = ncol(l$scores))
      else
        res <- Matrix::chol(l$B)

      return(res)}))
  }

  if(verbose)
    cat("Calculating MFPCA (", format(Sys.time(), "%T"), ")\n", sep = "")
  
  mArgvals <- if (utils::packageVersion("funData") <= "1.2") {
    getArgvals(mFData)
  } else {
    funData::argvals(mFData)
  }

  res <- calcMFPCA(N = N, p = p, Bchol = Bchol, M = M, type = type, weights = weights,
                   npc = npc, argvals = mArgvals, uniBasis = uniBasis, fit = fit, approx.eigen = approx.eigen)
  
  res$meanFunction <- m # return mean function, too
  
  names(res$functions) <- names(mFData)
  
  if(fit)
  {
    res$fit <- m + res$fit # add mean function to fits
    names(res$fit) <- names(mFData)
  } 
  
  # give correct names
  namesList <- lapply(mFData, names)
  if(!all(vapply(namesList, FUN = is.null, FUN.VALUE = TRUE))) 
  {
    if(length(unique(namesList)) != 1)
      warning("Elements have different curve names. Use names of the first element for the results.")
    
    row.names(res$scores) <- namesList[[1]]
    
    if(fit)
      for(i in seq_len(p))
        names(res$fit[[i]]) <- namesList[[1]]
  }

  # bootstrap for eigenfunctions
  if(bootstrap)
  {
    if(verbose)
      cat("Bootstrapping results:\n")

    booteFuns <- vector("list", p)

    for(j in seq_len(p))
      booteFuns[[j]] <- array(NA, dim  = c(nBootstrap, M, vapply(mFData[[j]]@argvals, FUN = length, FUN.VALUE = 0)))
    
    booteVals <- matrix(NA, nrow = nBootstrap, ncol = M)

    for(n in seq_len(nBootstrap))
    {
      if(verbose)
      {
        if(n %% 10 == 0)
          cat("\t n = ", n, " (", format(Sys.time(), "%T"), ")\n", sep = "")
      }

      if(is.null(bootstrapStrat))
        bootObs <- sample(N, replace = TRUE)
      else
        bootObs <- stratSample(bootstrapStrat)

      bootBasis <- vector("list", p)

      for(j in seq_len(p))
      {
        if(!is.null(uniBasis[[j]]$functions)) # re-estimate scores AND functions
        {
          bootBasis[[j]] <- do.call(univDecomp, c(list(funDataObject = (mFData[bootObs])[[j]]), uniExpansions[[j]]))
          
          # recalculate Bchol if necessary
          if(!bootBasis[[j]]$ortho)
            Bchol[[j]] <- Matrix::chol(bootBasis[[j]]$B)
        } 
        else # resample scores (functions are given and scores can simply be resampled)
          bootBasis[[j]] <- list(scores = uniBasis[[j]]$scores[bootObs, ], B = uniBasis[[j]]$B, ortho = uniBasis[[j]]$ortho, functions = uniBasis[[j]]$functions,
                                 settings = uniBasis[[j]]$settings)
      }

      npcBoot <- vapply(bootBasis, function(x){dim(x$scores)[2]}, FUN.VALUE = 0) # get number of univariate basis functions

      if(M > sum(npcBoot))
        stop("Function MFPCA (bootstrap): total number of univariate basis functions must be greater or equal M!")

      # calculate MFPCA for bootstrap sample (Bchol has been updated for UMPCA)
      bootMFPCA <- calcMFPCA(N = N, p = p, Bchol = Bchol, M = M, type = type, weights = weights,
                           npc = npcBoot, argvals = mArgvals, uniBasis = bootBasis, fit = FALSE, approx.eigen = approx.eigen)

      # save eigenvalues
      booteVals[n,] <- bootMFPCA$values
      
      # flip bootstrap estimates if necessary
      tmpFuns <- flipFuns(res$functions, bootMFPCA$functions)

      # save in booteFuns
      for(j in seq_len(p))
      {
        if(dimSupp[j] == 1)
          booteFuns[[j]][n,,] <- tmpFuns[[j]]@X
        if(dimSupp[j] == 2)
          booteFuns[[j]][n,,,] <- tmpFuns[[j]]@X
        if(dimSupp[j] == 3)
          booteFuns[[j]][n,,,,] <- tmpFuns[[j]]@X
      }
    }

    CIvalues <- vector("list", length(bootstrapAlpha))
    CI <- vector("list", length(bootstrapAlpha))
   
    for(alpha in seq_len(length(bootstrapAlpha)))
    {
      if(verbose)
        cat("Calculating bootstrap quantiles for alpha = ", bootstrapAlpha[alpha], " (", format(Sys.time(), "%T"), ")\n", sep = "")

      CIvalues[[alpha]]$lower <- apply(booteVals, 2, quantile, bootstrapAlpha[alpha]/2)
      CIvalues[[alpha]]$upper <-  apply(booteVals, 2, quantile, 1 - bootstrapAlpha[alpha]/2)
      names(CIvalues)[alpha] <- paste("alpha", bootstrapAlpha[alpha], sep = "_")
      
      bootCI_lower <- bootCI_upper <-  vector("list", p)
      for(j in seq_len(p))
      {
        bootCI_lower[[j]] <- funData(mFData[[j]]@argvals, apply(booteFuns[[j]], 2:length(dim(booteFuns[[j]])),
                                                             quantile, bootstrapAlpha[alpha]/2))
        bootCI_upper[[j]] <- funData(mFData[[j]]@argvals, apply(booteFuns[[j]],  2:length(dim(booteFuns[[j]])),
                                                             quantile, 1 - bootstrapAlpha[alpha]/2))
      }

      CI[[alpha]]$lower <- multiFunData(bootCI_lower)
      CI[[alpha]]$upper <- multiFunData(bootCI_upper)

      names(CI)[alpha] <- paste("alpha", bootstrapAlpha[alpha], sep = "_")
    }
    
    res$CIvalues <- CIvalues

    res$CI <- CI
  }
  
  class(res) <- "MFPCAfit"

  return(res)
}

#' Internal function that implements the MFPCA algorithm for given univariate decompositions
#'
#' @importFrom stats cov
#' @importFrom Matrix t
#' @importFrom foreach foreach
#' @importFrom irlba irlba
#'
#' @keywords internal
calcMFPCA <- function(N, p, Bchol, M, type, weights, npc, argvals, uniBasis, fit = FALSE, approx.eigen = FALSE)
{
  # combine all scores
  allScores <- foreach::foreach(j = seq_len(p), .combine = "cbind")%do%{uniBasis[[j]]$scores}

  # block vector of weights
  allWeights <- foreach::foreach(j = seq_len(p), .combine = "c")%do%{rep(sqrt(weights[j]), npc[j])}

  Z <- allScores %*% Matrix::Diagonal(x = allWeights) / sqrt(N-1)
  
  # check if approximation is appropriate (cf. irlba)
  if(approx.eigen & (M > min(N, sum(npc))/2))
  {
    warning("Calculating a large percentage of principal components, approximation may not be appropriate.
            'approx.eigen' set to FALSE.")
    approx.eigen = FALSE
  }

  # check if non-orthonormal basis functions used and calculate PCA on scores
  if(is.null(Bchol))
  {
    if(approx.eigen)
    {
      tmpSVD <- irlba::irlba(as.matrix(Z), nv = M)

      vectors <- tmpSVD$v
      values <- tmpSVD$d[seq_len(M)]^2
    }
    else
    {
      if(sum(npc) > 1000)
        warning("MFPCA with > 1000 univariate eigenfunctions and approx.eigen = FALSE. This may take some time...")

      e <- eigen(stats::cov(allScores) * outer(allWeights, allWeights, "*"))

      values <- e$values[seq_len(M)]
      vectors <- e$vectors[,seq_len(M)]
    }
  }
  else
  {
    if(approx.eigen)
    {
      tmpSVD <- irlba::irlba(as.matrix(Matrix::tcrossprod(Z, Bchol)), nv = M)

      vectors <- Matrix::crossprod(Bchol, tmpSVD$v)
      values <- tmpSVD$d[seq_len(M)]^2
    }
    else
    {
      if(sum(npc) > 1000)
        warning("MFPCA with > 1000 univariate eigenfunctions and approx.eigen = FALSE. This may take some time...")

      e <- eigen(Matrix::crossprod(Bchol) %*% (stats::cov(allScores) * outer(allWeights, allWeights, "*")))

      values <- Re(e$values[seq_len(M)])
      vectors <- Re(e$vectors[,seq_len(M)])
    }
  }

  # normalization factors
  normFactors <- 1/sqrt(diag(as.matrix(Matrix::crossprod(Z %*% vectors))))

  ### Calculate scores
  scores <- Z %*% vectors * sqrt(N-1) # see defintion of Z above!
  scores <- as.matrix(scores %*% diag(sqrt(values) * normFactors, nrow = M, ncol = M)) # normalization

  ### Calculate eigenfunctions (incl. normalization)
  npcCum <- cumsum(c(0, npc)) # indices for blocks (-1)

  tmpWeights <- as.matrix(Matrix::crossprod(Z, Z %*%vectors))
  eFunctions <- foreach::foreach(j = seq_len(p)) %do% {
    univExpansion(type = type[j],
                  scores = 1/sqrt(weights[j] * values) * normFactors * t(tmpWeights[npcCum[j]+seq_len(npc[j]), , drop = FALSE]),
                  argvals = argvals[[j]],
                  functions = uniBasis[[j]]$functions,
                  params = uniBasis[[j]]$settings)
  }

  res <- list(values = values,
              functions = multiFunData(eFunctions),
              scores = scores,
              vectors = vectors,
              normFactors = normFactors)

  # calculate truncated Karhunen-Loeve representation (no mean here)
  if(fit)
    res$fit <- multivExpansion(multiFuns = res$functions, scores = scores)

  return(res)
}

