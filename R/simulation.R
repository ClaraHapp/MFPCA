#' @import funData
NULL

#' @importFrom methods setGeneric setMethod
NULL

#' Generate a sparse version of functional data objects
#'
#' This function generates an artificially sparsified version of a functional
#' data object of class \link[funData]{funData} (univariate) or
#' \link[funData]{multiFunData} (multivariate). The minimal and maximal number
#' of observation points for all observations can be supplied by the user.
#'
#' The technique for artificially sparsifying the data is as described in Yao et
#' al. (2005): For each element \eqn{x_i^{(j)}}{x_i^(j)} of an observed
#' (multivariate) functional data object \eqn{x_i}, a random number
#' \eqn{R_i^{(j)} \in \{\mathrm{minObs}, \ldots, \mathrm{maxObs}\}}{R_i^(j)
#' in {minObs, \ldots, maxObs}} of observation points is generated. The points
#' are sampled uniformly from the full grid \eqn{\{t_{j,1} , \ldots , t_{j,
#' S_j}\} \subset \mathcal{T}_j}{{t_{j,1} , \ldots , t_{j, S_j}} in \calT_j}, resulting in
#' observations \deqn{ x_{i,r}^{(j)} = x_i^{(j)}(t_{j,r}), \quad r = 1
#' ,\ldots,R_i^{(j)},~ j = 1, \ldots, p.}{ x_{i,r}^(j) = x_i^(j)(t_{j,r}), r = 1
#' ,\ldots,R_i^(j), j = 1, \ldots, p.}
#'
#' @param funDataObject A functional data object of class
#'   \code{\link[funData]{funData}} or \link[funData]{multiFunData}.
#' @param minObs The minimal number of observation points. Must be a scalar for
#'   univariate functional data (\code{\link[funData]{funData}} class) or a
#'   vector of the same length as \code{funDataObject} for multivariate
#'   functional data (\code{\link[funData]{multiFunData}} class), giving the
#'   minimal number of observations for each element. See Details.
#' @param maxObs The maximal number of observation points. Must be a scalar for
#'   univariate functional data (\code{\link[funData]{funData}} class) or a
#'   vector of the same length as \code{funDataObject} for multivariate
#'   functional data (\code{\link[funData]{multiFunData}} class), giving the
#'   minimal number of observations for each element. See Details.
#'
#' @return An object of the same class as \code{funDataObject}, which is a
#'   sparse version of the original data.
#'
#' @seealso \code{\link[funData]{funData}}, \code{\link[funData]{multiFunData}},
#'   \code{\link{simFunData}}, \code{\link{simMultiFunData}},
#'   \code{\link{addError}}.
#'
#' @references Yao, F., H.-G. Mueller and J.-L. Wang (2005) Functional Data
#'   Analysis for Sparse Longitudinal Data. Journal of the American Statistical
#'   Association, 100 (470), 577--590.
#'
#' @export sparsify
#'
#' @examples
#' oldPar <- par(no.readonly = TRUE)
#' set.seed(1)
#'
#' # univariate functional data
#' full <- simFunData(argvals = seq(0,1, 0.01), M = 10, eFunType = "Fourier",
#'                    eValType = "linear", N = 3)$simData
#' sparse <- sparsify(full, minObs = 4, maxObs = 10)
#'
#' plot(full, main = "Sparsify")
#' plot(sparse, type = "p", pch = 20, add = TRUE)
#' legend("topright", c("Full", "Sparse"), lty=  c(1, NA), pch = c(NA, 20))
#'
#' # Multivariate
#' full <- simMultiFunData(type = "split", argvals = list(seq(0,1, 0.01), seq(-.5,.5, 0.02)),
#'                         M = 10, eFunType = "Fourier", eValType = "linear", N = 3)$simData
#' sparse <- sparsify(full, minObs = c(4, 30), maxObs = c(10, 40))
#'
#' par(mfrow = c(1,2))
#' plot(full[[1]], main = "Sparsify (multivariate)", sub = "minObs = 4, maxObs = 10")
#' plot(sparse[[1]], type = "p", pch = 20, add = TRUE)
#'
#' plot(full[[2]], main = "Sparsify (multivariate)", sub= "minObs = 30, maxObs = 40")
#' plot(sparse[[2]], type = "p", pch = 20, add = TRUE)
#' legend("bottomright", c("Full", "Sparse"), lty=  c(1, NA), pch = c(NA, 20))
#'
#' par(oldPar)
setGeneric("sparsify", function(funDataObject, minObs, maxObs) {standardGeneric("sparsify")})


#' sparsify for univariate functional data
#'
#' @keywords internal
setMethod("sparsify", signature = "funData",
          function(funDataObject, minObs, maxObs){
            sparseData <- funDataObject

            n <- length(sparseData@argvals[[1]]) # total number of observation points

            for(i in 1:nObs(sparseData)) # for all observed functions
            {
              Ni <- sample(minObs:maxObs, 1) # number of observation points
              notNA <- sample(n, Ni) # sample observation points
              sparseData@X[i, -notNA] <- NA  # set all other values to NA
            }
            return(sparseData)})


#' sparsify for multivariate functional data
#'
#' @keywords internal
setMethod("sparsify", signature = "multiFunData",
          function(funDataObject, minObs, maxObs){
            return(multiFunData(mapply(function(dat, minObs, maxObs){sparsify(dat, minObs, maxObs)}, funDataObject, minObs, maxObs)))
          })


#' Add Gaussian white noise to functional data objects
#'
#' This function generates an artificial noisy version of a functional data
#' object of class \link[funData]{funData} (univariate) or
#' \link[funData]{multiFunData} (multivariate) by adding iid. realizations of
#' Gaussian random variables \eqn{\varepsilon \sim \mathcal{N}(0, \sigma^2)}{\eps ~ N(0, \sigma^2)} to
#' the observations. The standard deviation \eqn{\sigma} can be supplied by the
#' user.
#'
#' @param funDataObject A functional data object of class
#'   \code{\link[funData]{funData}} or \link[funData]{multiFunData}.
#' @param sd The standard deviation \eqn{\sigma} of the Gaussian white noise
#'   that is added to the data. Defaults to 1. See Description.
#'
#' @return An object of the same class as \code{funDataObject}, which is a noisy
#'   version of the original data.
#'
#' @seealso \code{\link[funData]{funData}},
#'   \code{\link[funData]{multiFunData}}, \code{\link{simFunData}},
#'   \code{\link{simMultiFunData}}.
#'
#' @export addError
#'
#' @examples
#' oldPar <- par(no.readonly = TRUE)
#' set.seed(1)
#'
#' # Univariate functional data
#' plain <- simFunData(argvals = seq(0,1, 0.01), M = 10, eFunType = "Fourier",
#'                     eValType = "linear", N = 1)$simData
#' noisy <- addError(plain , sd = 0.5)
#' veryNoisy <- addError(plain, sd = 2)
#'
#' plot(plain, main = "Add error", ylim = range(veryNoisy@@X))
#' plot(noisy, type = "p", pch = 20, add = TRUE)
#' plot(veryNoisy, type = "p", pch = 4, add = TRUE)
#' legend("topright", c("Plain", "Noisy", "Very Noisy"), lty=  c(1, NA, NA), pch = c(NA, 20 ,4))
#'
#' # Multivariate functional data
#' plain <- simMultiFunData(type = "split", argvals = list(seq(0,1, 0.01), seq(-.5,.5, 0.02)), M = 10,
#'                         eFunType = "Fourier", eValType = "linear", N = 1)$simData
#' noisy <- addError(plain , sd = 0.5)
#' veryNoisy <- addError(plain, sd = 2)
#'
#' par(mfrow = c(1,2))
#' plot(plain[[1]], main = "Add error (multivariate)", ylim = range(veryNoisy[[1]]@@X))
#' plot(noisy[[1]], type = "p", pch = 20, add = TRUE)
#' plot(veryNoisy[[1]], type = "p", pch = 4, add = TRUE)
#'
#' plot(plain[[2]], main = "Add error (multivariate)", ylim = range(veryNoisy[[2]]@@X))
#' plot(noisy[[2]], type = "p", pch = 20, add = TRUE)
#' plot(veryNoisy[[2]], type = "p", pch = 4, add = TRUE)
#' legend("topright", c("Plain", "Noisy", "Very Noisy"), lty=  c(1, NA, NA), pch = c(NA, 20 ,4))
#'
#' par(oldPar)
setGeneric("addError", function(funDataObject, sd) {standardGeneric("addError")})


#'  Add gaussian white noise to functional data
#'
#' @keywords internal
setMethod("addError", signature = "funData",
          function(funDataObject, sd){
            ME <- array(rnorm(prod(dim(funDataObject@X)), mean = 0, sd = sd),  dim(funDataObject@X))

            return(funDataObject + funData(funDataObject@argvals, ME))
          })

#' Add gaussian white noise to multivariate functional data
#'
#' @keywords internal
setMethod("addError", signature = "multiFunData",
          function(funDataObject, sd){
            return(multiFunData(mapply(function(dat, sd){addError(dat, sd)}, funDataObject, sd)))
          })


#' Legendre Polynomials of degree 0,...,M-1
#'
#' This function iteratively calculates orthonormal Legendre polynomials of
#' degree 0,...,M-1 on an arbitrary interval.
#'
#' @param argvals A vector, defining a (fine) grid on the interval for which the
#'   Legendre polynomials are computed.
#' @param M An integer, specifying the number (and hence the degree) of
#'   polynomials that are calculcated.
#'
#' @return A univariate functional data object of class \code{\link[funData]{funData}}
#'   containing the Legendre polynomials on the given interval.
#'
#' @seealso \code{\link[funData]{funData}}, \code{\link{simFunData}}, \code{\link{simMultiFunData}}.
#'
#' @keywords internal
efPoly <- function(argvals, M)
{
  Phi <- matrix(NA, ncol = length(argvals), nrow = M)

  Phi[1, ] <- 1
  if(M == 1)
    return(funData(argvals, Phi))

  Phi[2, ] <- (2*(argvals-min(argvals))/diff(range(argvals))-1)

  if(M > 2)
  {
    for(m in 3:M) # for each function
      Phi[m, ] <- (2*m-3)/(m-1)*(2*(argvals-min(argvals))/diff(range(argvals))-1)*Phi[m-1, ]- (m-2)/(m-1)*Phi[m-2, ]
  }

  for(m in 1:M) # normalize
    Phi[m, ] <- sqrt((2*m-1)/diff(range(argvals)))*Phi[m, ]

  return(funData(argvals, Phi))
}


#' Calculate the first M Fourier basis functions
#'
#' This function calculates the first M orthonormal Fourier basis functions on
#' an arbitrary interval.
#'
#' If \code{linear}, the last basis function does not belong to the  Fourier
#' basis, but is the linear function orthogonalized to all previous Fourier
#' basis functions via the Gram-Schmidt method. This is implemented only if
#' \code{argvalss} is a grid defining the unit interval \eqn{[0,1]}.
#'
#' @param argvals A vector, defining a (fine) grid on the interval for which the
#'   Fourier basis functions are computed.
#' @param M An integer, specifying the number of basis functions that are
#'   calculcated.
#' @param linear Logical. If \code{TRUE}, the last function is not a Fourier
#'   function but the linear function orthogonalized to all previous Fourier
#'   basis functions. Defaults to \code{FALSE}. See Details.
#'
#' @return A univariate functional data object of class
#'   \code{\link[funData]{funData}} containing the Fourier basis functions on
#'   the given interval.
#'
#' @seealso \code{\link[funData]{funData}}, \code{\link{simFunData}}, \code{\link{simMultiFunData}}.
#'
#' @keywords internal
efFourier <- function(argvals, M, linear=F)
{
  Phi <- matrix(NA, nrow = M, ncol = length(argvals))

  Phi[1,] <- sqrt(1/diff(range(argvals)))

  if(M == 1)
    return(funData(argvals, Phi))

  for(m in 2:M)
  {
    if(m%%2 == 0) # m even
      Phi[m, ] <- sqrt(2/diff(range(argvals)))*cos((m%/%2)*(2*pi*(argvals-min(argvals))/diff(range(argvals))-pi))
    else # m odd
      Phi[m, ] <- sqrt(2/diff(range(argvals)))*sin((m%/%2)*(2*pi*(argvals-min(argvals))/diff(range(argvals))-pi))
  }

  if(linear) # overwrite Phi[M, ], add linear function and orthonormalize (Gram-Schmidt)
  {
    if(!all.equal(range(argvals) , c(0,1)))
      stop("efFourier, option linear: not yet implemented for argvalss != [0,1]!")

    # orthogonalize (exact function)
    Phi[M, ] <- argvals - 1/2  + rowSums(apply(matrix(1:((M-1)%/%2)), 1, function(k) (-1)^k/(pi*k)*sin(k*(2*pi*argvals-pi))))
    # normalize
    Phi[M, ] <-  Phi[M, ]/sqrt(1/3-1/4 -1/(2*pi^2)* sum( 1/(1:((M-1)%/%2))^2 ))
  }

  return(funData(argvals, Phi))
}


#' Calculate the first M eigenfunctions of the Wiener process
#'
#' This function calculates the first M (orthonormal) eigenfunctions of the
#' Wiener process on an arbitrary interval.
#'
#' @param argvals A vector, defining a (fine) grid on the interval for which the
#'   eigenfunctions are computed.
#' @param M An integer, specifying the number of eigenfunctions that are
#'   calculcated.
#'
#' @return A univariate functional data object of class
#'   \code{\link[funData]{funData}} containing the eigenfunctions of the Wiener
#'   process on the given interval.
#'
#' @seealso \code{\link[funData]{funData}}, \code{\link{simFunData}}, \code{\link{simMultiFunData}}.
#'
#' @keywords internal
efWiener <- function(argvals, M)
{
  Phi <- sapply(1:M, function(m,t){sqrt(2/diff(range(t)))*sin( (pi/2) * (2*m - 1) * (t-min(t))/diff(range(t)))}, t=argvals)

  return(funData(argvals, t(Phi)))
}


#' Generate orthonormal eigenfunctions
#'
#' This function calculates \eqn{M} (orthonormal) basis functions on a given
#' interval, that can be interpreted as the first \eqn{M} eigenfunctions of an
#' appropriate data generating process of functional data.
#'
#' The function implements three families of orthonormal basis functions plus
#' variations of them. The parameter \code{type}, that specifies the functions
#' to be calculated, can have the following values: \itemize{\item \code{"Poly"}:
#' Calculate orthonormal Legendre polynomials of degree 0,...,M-1. \item
#' \code{"PolyHigh"}: Calculate \eqn{M} orthonormal Legendre Polynomials of
#' higher degree. The vector of indices \code{ignoreDeg} specifies the functions
#' to be ignored. If \code{ignoreDeg} is not specified, the function returns an
#' error. \item \code{"Fourier"}: Calculate the first \eqn{M} Fourier basis functions.
#' \item \code{"FourierLin"}: Calculate the first \eqn{M-1} Fourier basis functions
#' plus the linear function, orthonormalized to the previous funcitons via
#' Gram-Schmidts method. This type is currently implemented for functions on the
#' unit interval \eqn{[0,1]} only. \item \code{"Wiener"}: Calculate the first \eqn{M}
#' orthonormal eigenfunctions of the Wiener process. }
#'
#' @param argvals A vector of numerics, defining a (fine) grid on the interval for which the
#'   basis functions are computed.
#' @param M An integer, specifying the number of functions that are
#'   calculated.
#' @param ignoreDeg A vector of numerics, specifying the degrees to be ignored
#'   for type \code{"PolyHigh"}. Defaults to \code{NULL}. See Details.
#' @param type A character string, specifying the type of functions that
#'   are calculated. See Details.
#'
#' @return A univariate functional data object of class
#'   \code{\link[funData]{funData}} containing the basis functions on the given
#'   interval.
#'
#' @seealso \link[funData]{funData}, \code{\link{simFunData}}, \code{\link{simMultiFunData}}.
#'
#' @export eFun
#'
#'  @examples
#' oldPar <- par(no.readonly = TRUE)
#'
#' argvals <- seq(0,1,0.01)
#'
#' par(mfrow = c(3,2))
#' plot(eFun(argvals, M = 4, type = "Poly"), main = "Poly", ylim = c(-3,3))
#' plot(eFun(argvals, M = 4, ignoreDeg = 1:2, type = "PolyHigh"), main = "PolyHigh",  ylim = c(-3,3))
#' plot(eFun(argvals, M = 4, type = "Fourier"), main = "Fourier", ylim = c(-3,3))
#' plot(eFun(argvals, M = 4, type = "FourierLin"), main = "FourierLin", ylim = c(-3,3))
#' plot(eFun(argvals, M = 4, type = "Wiener"), main = "Wiener",  ylim = c(-3,3))
#' par(oldPar)
eFun <- function(argvals, M, ignoreDeg = NULL, type)
{
  ret <- switch(type,
                Poly = efPoly(argvals, M),
                PolyHigh ={
                  if(is.null(ignoreDeg ))
                    stop("eFun, type = PolyHigh: specify ignoreDeg !")

                  extractObs(efPoly(argvals, M+length(ignoreDeg)), obs = -ignoreDeg)
                },
                Fourier = efFourier(argvals, M, linear = FALSE),
                FourierLin = efFourier(argvals, M, linear = TRUE),
                Wiener = efWiener(argvals, M),
                stop("eFun: choose either Poly, PolyHigh, Fourier, FourierLin or Wiener"))
  return(ret)
}


#' Generate eigenvalues
#'
#' This function generates \eqn{M} decreasing eigenvalues.
#'
#' The function implements three types of eigenvalues: \itemize{\item
#' \code{"linear"} The eigenvalues start at  \eqn{1} and decrease linearly
#' towards \eqn{0}: \deqn{\nu_m = \frac{M+1 -m}{m}.}{\nu_m = (M+1-m)/m.} \item \code{"exponential"}
#' The eigenvalues start at \eqn{1}  and decrease exponentially towards \eqn{0}:
#' \deqn{\nu_m = \exp\left(-\frac{m-1}{2}\right).}{\nu_m = exp(-(m-1)/2).}\item \code{"wiener"} The eigenvalues
#' correspond to the eigenvalueso of the Wiener process: \deqn{\nu_m =
#' \frac{1}{(\pi/2 \cdot (2m-1))^2}.}{\nu_m = (pi/2 * (2m-1))^(-2)} }
#'
#' @param M An integer, the number of eigenvalues to be generated.
#' @param type A character string specifying the type of eigenvalues that should
#'   be calculated. See Details.
#'
#' @return A vector containing the \eqn{M} decreasing eigenvalues.
#'
#' @export eVal
eVal <- function(M, type)
{
  ret <- switch(type,
                linear = ((M+1)-(1:M))/M,
                exponential = exp(-(0:(M-1))/2),
                wiener = 1/(pi/2 * (2*(1:M)-1))^2,
                stop("eVal: choose either linear, exponential or wiener"))
  return(ret)
}


#' Simulate univariate functional data
#'
#' This functions simulates (univariate) functional data \eqn{f_1,\ldots, f_N} based on a truncated
#' Karhunen-Lo\`{e}ve representation: \deqn{f_i(t) = \sum_{m = 1}^M \xi_{i,m}
#' \phi_m(t).} The eigenfunctions (basis functions) \eqn{\phi_m(t)} are generated
#' using \code{\link{eFun}}, the scores \eqn{\xi_{i,m}} are simulated
#' independently from a normal distribution with zero mean and decreasing
#' variance based on the \code{\link{eVal}} function.
#'
#' @param argvals A numeric vector, containing the observation points (a fine grid
#'   on a real interval) of the functional data that is to be simulated.
#' @param M An integer, giving the number of unvariate basis functions to use.
#'   See Details.
#' @param eFunType A character string specifying the type of univariate
#'   orthonormal basis functions to use. See \code{\link{eFun}} for details.
#' @param ignoreDeg A vector of integers, specifying the degrees to ignore when
#'   generating the univariate orthonormal bases. Defaults to \code{NULL}. See
#'   \code{\link{eFun}} for details.
#' @param eValType A character string, specifying the type of
#'   eigenvalues/variances used for the generation of the simulated functions
#'   based on the truncated Karhunen-Lo\`{e}ve representation. See
#'   \code{\link{eVal}} for details.
#' @param N An integer, specifying the number of multivariate functions to be
#'   generated.
#'
#' @return \item{simData}{A \code{\link[funData]{funData}} object with \code{N}
#'   observations, representing the simulated functional data.}
#'   \item{trueFuns}{A \code{\link[funData]{funData}} object with \code{M}
#'   observations, representing the true eigenfunction basis used for simulating
#'   the data.} \item{trueVals}{A vector of numerics, representing the true
#'   eigenvalues used for simulating the data.}
#'
#' @seealso \link[funData]{funData}, \code{\link{eFun}}, \code{\link{eVal}},
#'   \code{\link{addError}}, \code{\link{sparsify}}.
#'
#' @export simFunData
#'
#' @examples
#' oldPar <- par(no.readonly = TRUE)
#'
#' # Use Legendre polynomials as eigenfunctions and a linear eigenvalue decrease
#' test <- simFunData(seq(0,1,0.01), M = 10, eFunType = "Poly", eValType = "linear", N = 10)
#'
#' plot(test$trueFuns, main = "True Eigenfunctions")
#' plot(test$simData, main = "Simulated Data")
#'
#' # The use of ignoreDeg for eFunType = "PolyHigh"
#' test <- simFunData(seq(0,1,0.01), M = 4, eFunType = "Poly", eValType = "linear", N = 10)
#' test_noConst <-  simFunData(seq(0,1,0.01), M = 4, eFunType = "PolyHigh",
#'                             ignoreDeg = 1, eValType = "linear", N = 10)
#' test_noLinear <-  simFunData(seq(0,1,0.01), M = 4, eFunType = "PolyHigh",
#'                              ignoreDeg = 2, eValType = "linear", N = 10)
#' test_noBoth <-  simFunData(seq(0,1,0.01), M = 4, eFunType = "PolyHigh",
#'                            ignoreDeg = 1:2, eValType = "linear", N = 10)
#'
#' par(mfrow = c(2,2))
#' plot(test$trueFuns, main = "Standard polynomial basis (M = 4)")
#' plot(test_noConst$trueFuns, main = "No constant basis function")
#' plot(test_noLinear$trueFuns, main = "No linear basis function")
#' plot(test_noBoth$trueFuns, main = "Neither linear nor constant basis function")
#'
#' par(oldPar)
simFunData <- function(argvals, M, eFunType, ignoreDeg = NULL, eValType, N)
{
  # generate eigenfunctions
  trueFuns <-  eFun(argvals, M, ignoreDeg = ignoreDeg, type = eFunType)

  # generate eigenvalues and scores
  trueVals <- eVal(M, eValType)
  scores <- t(replicate(N, rnorm(M, sd = sqrt(eVal(M, eValType)))))

  # truncated Karhunen-Lo\`{e}ve representation
  simData <- funData(argvals, scores %*% trueFuns@X)

  return(list(simData = simData,
              trueFuns = trueFuns,
              trueVals = trueVals))

}


#' Simulate multivariate functional data
#'
#' This function provides a unified simulation structure for multivariate
#' functional data \eqn{f_1, \ldots, f_N} on one- or two-dimensional domains,
#' based on a truncated multivariate Karhunen-Lo\`{e}ve representation:
#' \deqn{f_i(t) = \sum_{m = 1}^M \rho_{i,m} \psi_m(t).} The multivariate
#' eigenfunctions (basis functions) \eqn{\psi_m}  are constructed from
#' univariate orthonormal bases. There are two different concepts for the
#' construction, that can be chosen by the parameter \code{type}: A split
#' orthonormal basis (\code{split}, only one-dimensional domains) and weighted
#' univariate orthonormal bases (\code{weighted}, one- and two-dimensional
#' domains). The scores \eqn{\rho_{i,m}} in the Karhunen-Lo\`{e}ve
#' representation are simulated independently from a normal distribution with
#' zero mean and decreasing variance. See Details.
#'
#' The parameter \code{type} defines how the eigenfunction basis for the
#' multivariate Karhunen-Lo\`{e}ve representation is constructed: \itemize{
#' \item \code{type = "split"}: The basis functions of an underlying 'big'
#' orthonormal basis are split in \code{M} parts, translated and possibly
#' reflected. This yields an orthornormal basis of multivariate functions with
#' \code{M} elements. This option is implemented only for one-dimensional
#' domains. \item \code{type = "weighted"} The multivariate eigenfunction basis
#' consists of weighted univariate orthonormal bases.  This yields an
#' orthornormal basis of multivariate functions with \code{M} elements. For data
#' on two-dimensional domains (images), the univariate basis is constructed as a
#' tensor product of univariate bases in each direction (x- and y-direction). }
#'
#' Depending on \code{type}, the other parameters have to be specified as
#' follows: \subsection{Split 'big' orthonormal basis:}{ The parameters \code{M}
#' (integer), \code{eFunType} (character string) and \code{ignoreDeg} (integer
#' vector or \code{NULL}) are passed to the function \code{\link{eFun}} to
#' generate a univariate orthornormal basis on a 'big' interval. Subsequently,
#' the basis functions are split and translated, such that the \eqn{j}-th part
#' of the split function is defined on the interval corresponding to
#' \code{argvals[[j]]}. The elements of the multivariate basis functions are
#' given by these split parts of the original basis functions multiplied by a
#' random sign \eqn{\sigma_j \in \{-1,1\}}{\sigma_j in {-1,1}}.}
#'
#' \subsection{Weighted orthonormal bases:}{ The parameters \code{argvals, M,
#' eFunType} and \code{ignoreDeg} are all lists of a similar structure. They are
#' passed element-wise to the function \code{\link{eFun}} to generate
#' orthornormal basis functions for each element of the multivariate functional
#' data to be simulated. In case of bivariate elements (images), the
#' corresponding basis functions are constructed as tensor products of
#' orthonormal basis functions in each direction (x- and y-direction).
#'
#' If the \eqn{j}-th element of the simulated data should be defined on a
#' one-dimensional domain, then \itemize{ \item \code{argvals[[j]]} is a list,
#' containing one vector of observation points. \item \code{M[[j]]} is an
#' integer, specifiying the number of basis functions to use for this entry.
#' \item  \code{eFunType[[j]]} is a character string, specifying the type of
#' orthonormal basis functions to use for this entry (see \code{\link{eFun}} for
#' possible options). \item \code{ignoreDeg[[j]]} is a vector of integers,
#' specifying the degrees to ignore when constructing the orthonormal basis
#' functions. The default value is \code{NULL}. }
#'
#' If the \eqn{j}-th element of the simulated data should be defined on a
#' two-dimensional domain, then \itemize{ \item \code{argvals[[j]]} is a list,
#' containing two vector of observation points, one for each direction
#' (observation points in x-direction and in y-direction). \item \code{M[[j]]}
#' is a vector of two integers, giving the number of basis functions for each
#' direction (observation points in x-direction and in y-direction). \item
#' \code{eFunType[[j]]} is a vector of two character strings, giving the type of
#' orthonormal basis functions for each direction (x- and y-direction, see
#' \link{eFun} for possible options). The corresponding basis functions are
#' constructed as tensor products of orthonormal basis functions in each
#' direction. \item \code{ignoreDeg[[j]]} is a list, containing two integer
#' vectors that specify the degrees to ignore when constructing the orthonormal
#' basis functions in each direction. The default value is \code{NULL}. } The
#' total number of basis functions (i.e. the product of \code{M[[j]]} for all
#' \code{j}) must be equal!}
#'
#' @param type A character string, specifying the construction method for the
#'   multivariate eigenfunctions (either \code{"split"} or \code{"weighted"}).
#'   See Details.
#' @param argvals A list, containing the observation points for each element of
#'   the multivariate functional data that is to be simulated. The length of
#'   \code{argvals} determines the number of elements in the resulting simulated
#'   multivariate functional data. See Details.
#' @param M An integer (\code{type = "split"}) or a list of integers (\code{type
#'   = "weighted"}), giving the number of unvariate basis functions to use. See
#'   Details.
#' @param eFunType A character string (\code{type = "split"})   or a list of
#'   character strings(\code{type = "weighted"}), specifying the type of
#'   univariate orthonormal basis functions to use. See Details.
#' @param ignoreDeg A vector of integers (\code{type = "split"})   or a list of
#'   integer vectors (\code{type = "weighted"}), specifying the degrees to
#'   ignore when generating the univariate orthonormal bases. Defaults to
#'   \code{NULL}. See Details.
#' @param eValType A character string, specifying the type of
#'   eigenvalues/variances used for the simulation of the multivariate functions
#'   based on the truncated Karhunen-Lo\`{e}ve representation. See
#'   \code{\link{eVal}} for details.
#' @param N An integer, specifying the number of multivariate functions to be
#'   generated.
#'
#' @return \item{simData}{A \code{\link[funData]{multiFunData}} object with
#'   \code{N} observations, representing the simulated multivariate functional
#'   data.} \item{trueFuns}{A \code{\link[funData]{multiFunData}} object with
#'   \code{M} observations, representing the multivariate eigenfunction basis
#'   used for simulating the data.} \item{trueVals}{A vector of numerics,
#'   representing the eigenvalues used for simulating the data.}
#'
#' @seealso \code{\link[funData]{multiFunData}}, \code{\link{eFun}},
#'   \code{\link{eVal}}, \code{\link{simFunData}}, \code{\link{addError}},
#'   \code{\link{sparsify}}.
#'
#' @references C. Happ, S. Greven (2015): Multivariate Functional Principal
#'   Component Analysis for Data Observed on Different (Dimensional) Domains.
#'   Preprint on arXiv: \url{http://arxiv.org/abs/1509.02029}
#'
#' @export simMultiFunData
#'
#' @examples
#' oldPar <- par(no.readonly = TRUE)
#'
#' # split
#' split <- simMultiFunData(type = "split", argvals = list(seq(0,1,0.01), seq(-0.5,0.5,0.02)),
#'                  M = 5, eFunType = "Poly", eValType = "linear", N = 7)
#'
#' par(mfrow = c(1,2))
#' plot(split$trueFuns, main = "Split: True Eigenfunctions", ylim = c(-2,2))
#' plot(split$simData, main = "Split: Simulated Data")
#'
#' # weighted (one-dimensional domains)
#' weighted1D <- simMultiFunData(type = "weighted",
#'                  argvals = list(list(seq(0,1,0.01)), list(seq(-0.5,0.5,0.02))),
#'                  M = c(5,5), eFunType = c("Poly", "Fourier"), eValType = "linear", N = 7)
#'
#' plot(weighted1D$trueFuns, main = "Weighted (1D): True Eigenfunctions", ylim = c(-2,2))
#' plot(weighted1D$simData, main = "Weighted (1D): Simulated Data")
#'
#' # weighted (one- and two-dimensional domains)
#' weighted <- simMultiFunData(type = "weighted",
#'                argvals = list(list(seq(0,1,0.01), seq(0,10,0.1)), list(seq(-0.5,0.5,0.01))),
#'                M = list(c(5,4), 20), eFunType = list(c("Poly", "Fourier"), "Wiener"),
#'                eValType = "linear", N = 7)
#'
#' plot(weighted$trueFuns, main = "Weighted: True Eigenfunctions (m = 2)", obs = 2)
#' plot(weighted$trueFuns, main = "Weighted: True Eigenfunctions (m = 15)", obs = 15)
#' plot(weighted$simData, main = "Weighted: Simulated Data (1st observation)", obs = 1)
#' plot(weighted$simData, main = "Weighted: Simulated Data (2nd observation)", obs = 2)
#'
#' par(oldPar)
simMultiFunData <- function(type, argvals, M, eFunType, ignoreDeg = NULL, eValType, N)
{
  # generate eigenfunctions
  trueFuns <- switch(type,
                split = simMultiSplit(argvals, M, eFunType, ignoreDeg, eValType, N),
                weighted = simMultiWeight(argvals, M, eFunType, ignoreDeg, eValType, N),
                stop("simMultiFunData: choose either 'split' or 'weighted' for the simulation of multivariate functional data.")
                )

  # number of eigenfunctions generated
  Mtotal <- nObs(trueFuns)

  # number of elements in multivariate functional basis
  p <- length(trueFuns)

  # generate eigenvalues and scores
  trueVals <- eVal(Mtotal, eValType)
  scores <- t(replicate(N, rnorm(Mtotal, sd = sqrt(eVal(Mtotal, eValType)))))

  # generate individual observations
  simData  <- vector("list", p)

  for(j in 1:p)
  {
    simData[[j]] <- funData(trueFuns[[j]]@argvals, array(0, dim = c(N, dim(trueFuns[[j]]@X)[-1])))

    if(dimSupp(trueFuns[[j]]) == 1)
      simData[[j]]@X <- scores %*% trueFuns[[j]]@X
    else # dimSupp(trueFuns[[j]]) == 2
    {
      for(i in 1:N)
      {
        for(m in 1:Mtotal)
          simData[[j]]@X[i,,] <- simData[[j]]@X[i,,] + scores[i,m]* trueFuns[[j]]@X[m,,]
      }
    }
  }

  return(list(simData = multiFunData(simData),
              trueFuns = trueFuns,
              trueVals = trueVals))
}


#' Simulate multivariate eigenfunctions based on a split 'big' ONB
#'
#' @keywords internal
simMultiSplit <- function(argvals, M, eFunType, ignoreDeg = NULL, eValType, N)
{
  # consistency check
  if(any( c(length(M), length(eFunType), length(eValType)) != 1) )
    stop("simMultiSplit: argvals, M, eFunType, eValType must all be of length 1!")

  # number of elements
  p <- length(argvals)

  # "rearrange" argvalss
  x <- vector("list", length = length(argvals))
  splitVals <- rep(NA, length(argvals)+1)

  x[[1]] <- unlist(argvals[[1]]) # convert to vector, if argvals[[1]] is a list
  splitVals[1:2] <- c(0, length(x[[1]]))

  for(i in 2:p)
  {
    x[[i]] <- unlist(argvals[[i]]) # convert to vector, if argvals[[i]] is a list
    x[[i]] <- argvals[[i]] - min(argvals[[i]]) + max(x[[i-1]])
    splitVals[i+1] <- splitVals[i]+length(x[[i]])
  }

  # generate "big" orthonormal system
  f <-  eFun(unlist(x), M, ignoreDeg = ignoreDeg, type = eFunType)

  # sample sign randomly
  s <- sample(c(-1,1), p, 0.5)

  # result object
  trueFuns  <- vector("list", p)

  for(j in 1:p)
    trueFuns[[j]] <- funData(argvals[[j]],  s[j] * f@X[,(1+splitVals[j]):splitVals[j+1]])

  return(multiFunData(trueFuns))
}




#' Simulate multivariate eigenfunctions based on weighted orthonormal bases
#'
#'@importFrom foreach %do%
#'
#' @keywords internal
simMultiWeight <- function(argvals, M, eFunType, ignoreDeg = NULL, eValType, N)
{
  p <- length(argvals)

  # dimension for each component
  dimsSupp <- foreach::foreach(j = 1:p, .combine = "c")%do%{length(argvals[[j]])}

  if(any(dimsSupp > 2))
    stop("Function simMultiWeight: method is not implemented for objects of dimension > 2!")

  if(p > 1)
  {
    if(do.call(all.equal, lapply(M, prod)))
    {
      Mtotal <- prod(M[[1]])
    }
    else
      stop("Function simMultiWeight: basis dimensions must be equal!")
  }
  else
  {
    Mtotal <- prod(M[[1]])
  }

  # mixing parameters
  alpha <- runif(p, 0.2, 0.8)
  weight <- sqrt(alpha/sum(alpha))

  # generate basis
  basis <- vector("list", p)

  for(j in 1:p)
  {
    if(dimsSupp[j] == 1) # one-dimensional
    {
      basis[[j]] <- weight[j]*eFun(argvals[[j]][[1]], M = M[[j]], ignoreDeg = ignoreDeg[[j]], type = eFunType[[j]])@X
    }
    else # dimsSupp[j] == 2, i.e. two-dimensional
    {
      # image basis
      y1 <- eFun(argvals[[j]][[1]], M = M[[j]][1], ignoreDeg = ignoreDeg[[j]][[1]], type = eFunType[[j]][1])@X
      y2 <- eFun(argvals[[j]][[2]], M = M[[j]][2], ignoreDeg = ignoreDeg[[j]][[2]], type = eFunType[[j]][2])@X

      basis[[j]] <- array(0, c(prod(M[[j]]), length(argvals[[j]][[1]]), length(argvals[[j]][[2]])))

      for(l in 1:M[[j]][1])
      {
        for(k in 1:M[[j]][2])
          basis[[j]][(l-1)*M[[j]][2]+k,,] <- weight[j]*y1[l,] %o% y2[k,]
      }
    }
  }

  trueFuns <- vector("list", p)

  for(j in 1:p)
    trueFuns[[j]] <- funData(argvals[[j]], basis[[j]])

  return(multiFunData(trueFuns))
}
