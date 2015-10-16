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
univDecomp <- function(type, data, params)
{
  params$funDataObject <- data
  
  res <- switch(type,
                "uFPCA" = do.call(),
                "splines1D" = do.call(splineBasis1D, params),
                "splines1Dpen" = do.call(splineBasis1Dpen, params),
                "splines2D" = do.call(splineBasis2D, params),
                "splines2Dpen" = do.call(splineBasis2Dpen, params),
                "DCT2D" = ...,
                stop("Univariate Decomposition for 'type' = ", type, " not defined!")
  )
  
  if(res$ortho == FALSE & is.null(res$B))
    stop("UnivDecomp: must provide integral matrix B for non-orthonormal basis functions.")
  
  return(list(scores = res$scores,
              B = ifelse(res$ortho == TRUE, NULL, res$B),
              ortho = res$ortho,
              functions = ifelse(type == "uFPCA", res$functions, NULL)
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
#'   containing the observed functional data samples and for which the bases 
#'   decomposition is calculated.
#' @param bs A character string, specifying the type of basis functions to be
#'   used. Please refer to \code{\link[mgcv]{smooth.terms}} for a list of
#'   possible basis functions.
#' @param m A numeric, the order of the spline basis. See  \code{\link[mgcv]{s}}
#'   for details.
#' @param k A numeric, the number of basis functions used.  See 
#'   \code{\link[mgcv]{s}} for details.
#'   
#' @return \item{scores}{A matrix of weights with dimension \code{N x k},
#'   reflecting the weights for each basis function in each observation.} 
#'   \item{B}{A matrix containing the scalar product of all pairs of basis
#'   functions.} \item{ortho}{Logical, set to \code{FALSE}, as basis functions
#'   are not orthogonal} \item{functions}{\code{NULL}, as basis functions are
#'   known}
#'   
#' @seealso univDecomp
#'   
#' @importFrom mgcv gam
splineBasis1D <- function(funDataObject, bs, m, k)
{
  N <- nObs(funDataObject)
  
  x <- funDataObject@xVal[[1]]
  
  # spline design matrix via gam
  desMat <- mgcv::gam(funDataObject@X[1, ] ~ s(x, bs = bs, m = m, k = k), fit = FALSE)$X
  
  # weights via lm -> no penalization
  scores <- t(apply(funDataObject@X, 1, function(f, dM){lm(f ~ dM - 1)$coef}, dM = desMat))
  
  return(list(scores = scores,
              B = .calcBasisIntegrals(t(desMat), 1, funDataObject@xVal),
              ortho = FALSE,
              functions = NULL
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
#'   containing the observed functional data samples and for which the bases 
#'   representation is calculated.
#' @param bs A character string, specifying the type of basis functions to be
#'   used. Please refer to \code{\link[mgcv]{smooth.terms}} for a list of
#'   possible basis functions.
#' @param m A numeric, the order of the spline basis. See  \code{\link[mgcv]{s}}
#'   for details.
#' @param k A numeric, the number of basis functions used.  See 
#'   \code{\link[mgcv]{s}} for details.
#' @param parallel Logical. If \code{TRUE}, the coefficients for the basis 
#'   functions are calculated in parallel. The implementation is based on the 
#'   \code{\link[foreach]{foreach}} function and requires a parallel backend
#'   that must be registered before. See \code{\link[foreach]{foreach}} for
#'   details.
#'   
#' @return \item{scores}{A matrix of weights with dimension \code{N x k}, 
#'   reflecting the weights for each basis function in each observation.} 
#'   \item{B}{A matrix containing the scalar product of all pairs of basis 
#'   functions.} \item{ortho}{Logical, set to \code{FALSE}, as basis functions 
#'   are not orthogonal} \item{functions}{\code{NULL}, as basis functions are 
#'   known}
#'   
#' @seealso univDecomp
#'   
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#' @importFrom mgcv gam
splineBasis1Dpen <- function(funDataObject, bs, m, k, parallel = FALSE)
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
  
  scores <- rbind(scores, g$coef)
  
  return(list(scores = scores,
              B = .calcBasisIntegrals(t(model.matrix(g)), 1, funDataObject@xVal),
              ortho = FALSE,
              functions = NULL
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
#'   containing the observed functional data samples and for which the bases 
#'   representation is calculated.
#' @param bs A vector of character strings (or a single character string),
#'   specifying the type of basis functions to be used. Please refer to
#'   \code{\link[mgcv]{te}} for a list of possible basis functions.
#' @param m A numeric vector (or a single number), the order of the spline
#'   basis. See \code{\link[mgcv]{s}} for details.
#' @param k An numeric vector (or a single number), the number of basis
#'   functions used.  See  \code{\link[mgcv]{s}} for details.
#'   
#' @return \item{scores}{A matrix of weights with dimension \code{N x K}, 
#'   reflecting the weights for each basis function in each observation, where 
#'   \code{K} is the total number of basis functions used.} \item{B}{A matrix 
#'   containing the scalar product of all pairs of basis functions.} 
#'   \item{ortho}{Logical, set to \code{FALSE}, as basis functions are not 
#'   orthogonal.} \item{functions}{\code{NULL}, as basis functions are known.}
#'   
#' @seealso univDecomp
#'   
#' @importFrom foreach %do%
#' @importFrom mgcv gam
splineBasis2D <- function(funDataObject, bs, m, k)
{
  N <- nObs(funDataObject)
  
  coord <- expand.grid(x = funDataObject@xVal[[1]], y = funDataObject@xVal[[2]])
  
  browser()
  
  # spline design matrix via gam
  desMat <- mgcv::gam(as.vector(funDataObject@X[1,,]) ~ te(coord$x, coord$y, bs = bs, m = m, k = k), data = coord, fit = FALSE)$X
  
  # weights via lm -> no penalization
  scores <- t(apply(funDataObject@X, 1, function(f, dM){lm(as.vector(f) ~ dM - 1)$coef}, dM = desMat))
  
  # extract basis functions (in the correct dimensions)
  B <- aperm(array(desMat, c(nObsPoints(funDataObject), ncol(scores))), c(3,1,2))

  return(list(scores = scores,
              B = .calcBasisIntegrals(B, 2, funDataObject@xVal),
              ortho = FALSE,
              functions = NULL
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
#'   containing the observed functional data samples and for which the bases 
#'   representation is calculated.
#' @param bs An vector of character strings (or a single character), the type of
#'   basis functions to be used. Please refer to \code{\link[mgcv]{te}} for a
#'   list of possible basis functions.
#' @param m A numeric vector (or a single number), the order of the spline
#'   basis. See \code{\link[mgcv]{s}} for details.
#' @param k A numeric vector (or a single number), the number of basis functions
#'   used.  See  \code{\link[mgcv]{s}} for details.
#' @param parallel Logical. If \code{TRUE}, the coefficients for the basis 
#'   functions are calculated in parallel. The implementation is based on the 
#'   \code{\link[foreach]{foreach}} function and requires a parallel backend 
#'   that must be registered before. See \code{\link[foreach]{foreach}} for 
#'   details.
#'   
#' @return \item{scores}{A matrix of weights with dimension \code{N x K}, 
#'   reflecting the weights for each basis function in each observation, where 
#'   \code{K} is the total number of basis functions used.} \item{B}{A matrix
#'   containing the scalar product of all pairs of basis functions.}
#'   \item{ortho}{Logical, set to \code{FALSE}, as basis functions are not
#'   orthogonal} \item{functions}{\code{NULL}, as basis functions are known}
#'   
#' @seealso univDecomp
#'   
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#' @importFrom mgcv bam
splineBasis2Dpen <- function(funDataObject, bs, m, k, parallel = FALSE)
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
  
  scores <- rbind(scores, g$coef)

  # extract basis functions (in the correct dimensions)
  B <- aperm(array(model.matrix(g), c(nObsPoints(funDataObject), ncol(scores))), c(3,1,2))

  return(list(scores = scores,
              B = .calcBasisIntegrals(B, 2, funDataObject@xVal),
              ortho = FALSE,
              functions = NULL
  ))
}