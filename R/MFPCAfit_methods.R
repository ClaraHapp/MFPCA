#' Scoreplot Generic
#'
#'  Redirects to \code{\link[graphics]{plot.default}}
#'
#' @param PCAobject A principal component object
#' @param ... Arguments passed from or to other methods
#' 
#' @export
scoreplot <- function(PCAobject, ...) UseMethod("scoreplot", PCAobject)
scoreplot.default <- function(PCAobject, ...) graphics::plot.default(PCAobject, ...)




#' Plot the Scores of a Multivariate Functional Principal Component Analysis
#' 
#' This function plots two scores of a multivariate functional principal 
#' component analysis for each observation.
#' 
#' @param PCAobject An object of class \code{MFPCAfit}, typically returned by the
#'   \link{MFPCA} function.
#' @param choices The indices of the scores that should by displayed. Defaults 
#'   to \code{1:2}, i.e. the scores corresponding to the two leading modes of 
#'   variability in the data.
#' @param scale Logical. Should the scores be scaled by the estimated
#'   eigenvalues to emphasize the proportions of total variance explained by the
#'   components. Defaults to \code{FALSE}.
#' @param ... Further parameters passed to the
#'   \code{\link[graphics]{plot.default}} function.
#'   
#' @return A bivariate plot of scores.
#' 
#' @seealso \code{\link{MFPCA}}
#' 
#' @export
#' 
#' @examples 
# Simulate multivariate functional data on one-dimensonal domains
#' # and calculate MFPCA (cf. MFPCA help)
#' set.seed(1)
#' # simulate data (one-dimensional domains)
#' sim <-  simMultiFunData(type = "split", argvals = list(seq(0,1,0.01), seq(-0.5,0.5,0.02)),
#'                        M = 5, eFunType = "Poly", eValType = "linear", N = 100)
#' # MFPCA based on univariate FPCA
#' PCA <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"),
#'                                                      list(type = "uFPCA")))
#' 
#' # Plot the first two scores
#' scoreplot(PCA) # no scaling (default)
#' scoreplot(PCA, scale = TRUE) # scale the scores by the first two eigenvalues 
scoreplot.MFPCAfit <- function(PCAobject, choices = 1:2, scale = FALSE, ...)
{
  if(!inherits(PCAobject,"MFPCAfit"))
    stop("Argument is not of class 'MFPCAfit'.")
  
  if(!(is.numeric(choices) & length(choices) == 2 & all(choices > 0)))
    stop("Parameter 'choices' must be a vector of length 2 with positive entries.")
  
  if(NCOL(PCAobject$scores) < max(choices))
    stop(paste("Argument choices requires", max(choices), 
               "scores, MFPCA object contains only", NCOL(PCAobject$scores), "."))
  
  if(!is.logical(scale))
    stop("Parameter 'scale' must be passed as a logical.")
  
  # check for labels, otherwise construct them TODO: rownames for scores in MFPCA
  lab <- {if(is.null(rownames(PCAobject$scores)))
    seq_len(NROW(PCAobject$scores))
    else
      rownames(PCAobject$scores)}
  
  # scale scores if scale = TRUE
  plotScore <- {
    if(scale) # multiply row-wise by eigenvalues
      t(PCAobject$values[choices] * t(PCAobject$scores[, choices]))
    else
      PCAobject$scores[, choices]
  }
  
  # plot scores
  graphics::plot.default(x = plotScore, y = NULL, type = "n", 
                         xlab = paste("Scores PC", choices[1]), ylab = paste("Scores PC", choices[2]), ...)
  graphics::text(x = plotScore, y = NULL, labels = lab, ...)
  
  graphics::abline(h = 0, lty = 3, col = "gray")
  graphics::abline(v = 0, lty = 3, col = "gray")
  
  invisible(NULL)
}

#' Screeplot for Multivariate Functional Principal Component Analysis
#'
#' This function plots the proportion of variance explained by the leading
#' eigenvalues in an MFPCA against the number of the principal component.
#'
#' @param x An object of class MFPCAfit, typically returned by a call to
#'   \code{\link{MFPCA}}.
#' @param npcs The number of eigenvalued to be plotted. Defaults to all
#'   eigenvalues if their number is less or equal to 10, otherwise show
#'   only the leading first 10 eigenvalues.
#' @param type The type of screeplot to be plotted. Can be either
#'   \code{"lines"} or \code{"barplot"}. Defaults to \code{"lines"}.
#' @param ylim The limits for the y axis. Can be passed either as a vector
#'   of length 2 or as \code{NULL} (default). In the second case,
#'   \code{ylim} is set to \code{(0,max(pve))}, with \code{pve} the
#'   proportion of variance explained by the principal
#'   components to be plotted.
#' @param main The title of the plot. Defaults to the variable name of
#'   \code{x}.
#' @param ... Other graphic parameters passed to
#'   \code{\link[graphics]{plot.default}} (for \code{type = "lines"}) or
#'   \code{\link[graphics]{barplot}} (for \code{type = "barplot"}).
#'
#' @seealso \code{\link{MFPCA}}, \code{\link[stats]{screeplot}}
#'
#' @importFrom stats screeplot
#' @importFrom graphics barplot
#' @importFrom graphics axis
#'
#' @export
#'
#' @examples
#' # Simulate multivariate functional data on one-dimensonal domains
#' # and calculate MFPCA (cf. MFPCA help)
#' set.seed(1)
#' # simulate data (one-dimensional domains)
#' sim <-  simMultiFunData(type = "split", argvals = list(seq(0,1,0.01), seq(-0.5,0.5,0.02)),
#'                        M = 5, eFunType = "Poly", eValType = "linear", N = 100)
#' # MFPCA based on univariate FPCA
#' PCA <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"),
#'                                                      list(type = "uFPCA")))
#'
#' # screeplot
#' screeplot(PCA) # default options
#' screeplot(PCA, npcs = 3, type = "barplot", main= "Screeplot")
screeplot.MFPCAfit <- function(x, npcs = min(10, length(x$values)), type = "lines", ylim = NULL,
                               main = deparse(substitute(x)), ...)
{
  if(!inherits(x, "MFPCAfit"))
    stop("Argument is not of class 'MFPCAfit'.")
  
  if(!all(is.numeric(npcs), length(npcs) == 1, npcs > 0, npcs <= length(x$values)))
    stop("Parameter 'npcs' must be a number between 1 and ", length(x$values), ".")
  
  if(!is.character(type))
    stop("Parameter 'type' must be passed as a character.")
  
  if(!(is.character(main) | is.null(main)))
    stop("Parameter 'main' must be either NULL or passed as a character.")

  ylab <- "Proportion of Variance Explained"
  pve <- x$values[seq_len(npcs)]/sum(x$values)
  
  if(is.null(ylim))
     ylim <- c(0,max(pve))
  
  switch(type,
         "lines" = {plot(x = seq_len(npcs), y = pve, type = "b", ylim = ylim,
                                   main = main, xlab = "PCs", ylab = ylab, xaxt = "n", ...)
                    graphics::axis(1, at = seq_len(npcs))},
         "barplot" = graphics::barplot(height = pve, main = main, ylim = ylim,
                                       names.arg = paste("PC", seq_len(npcs)), ylab = ylab, ...),
         stop("Type ", type, " not defined in screeplot.")
  )
  
  invisible(NULL)
}


#' Plot MFPCA results
#' 
#' Plots the eigenfunctions as perturbations of the mean (i.e. the mean function
#' plus/minus a constant factor times each eigenfunction separately). If all 
#' elements have a one-dimensional domain, the plots can be combined, otherwise 
#' the effects of adding and subtracting are shown in two separate rows for each
#' eigenfunction.
#' 
#' @param x An object of class \code{MFPCAfit}, typically returned by the 
#'   \link{MFPCA} function.
#' @param plotPCs The principal components to be plotted. Defaults to all 
#'   components in the \code{MFPCAfit} object.
#' @param stretchFactor The factor by which the principal components are 
#'   multiplied before adding / subtracting them from the mean function. If 
#'   \code{NULL} (the default), the median absolute value of the scores of each 
#'   eigenfunction is used.
#' @param combined Logical: Should the plots be combined? (Works only if all 
#'   dimensions are one-dimensional). Defaults to \code{FALSE}.
#' @param ... Further graphical parameters passed to the
#'   \link[funData]{plot.funData} functions for functional data.
#'      
#' @return A plot of the principal components as perturbations of the mean.
#'   
#' @seealso \code{\link{MFPCA}}, \code{\link[funData]{plot.funData}}
#'   
#' @method plot MFPCAfit 
#' @export
#' 
#' @examples 
#' # Simulate multivariate functional data on one-dimensonal domains
#' # and calculate MFPCA (cf. MFPCA help)
#' set.seed(1)
#' # simulate data (one-dimensional domains)
#' sim <-  simMultiFunData(type = "split", argvals = list(seq(0,1,0.01), seq(-0.5,0.5,0.02)),
#'                        M = 5, eFunType = "Poly", eValType = "linear", N = 100)
#' # MFPCA based on univariate FPCA
#' PCA <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"),
#'                                                      list(type = "uFPCA")))
#' 
#' # Plot the results
#' plot(PCA, combined = TRUE) # combine addition and subtraction in one plot
plot.MFPCAfit <- function(x, plotPCs = seq_len(nObs(x$functions)), stretchFactor = NULL, combined = FALSE, ...)
{
  if(!inherits(x, "MFPCAfit"))
    stop("Argument is not of class 'MFPCAfit'.")
  
  if(!(is.numeric(plotPCs) & 0 < length(plotPCs) & length(plotPCs) <= length(x$values) & all(0 < plotPCs, plotPCs <= length(x$values))))
    stop("Parameter 'plotPCs' must be a vector with values between 1 and ", length(x$values), ".")
  
  if(!(is.null(stretchFactor) | all(is.numeric(stretchFactor), length(stretchFactor) == 1, stretchFactor > 0)))
    stop("Parameter 'stretchFactor' must be either NULL or a positive number.")
  
  if(!is.logical(combined))
    stop("Parameter 'combined' must be passed as a logical.")
  
  # check dimensions
  dims <- funData::dimSupp(x$functions)
  
  if(any(dims > 2))
    stop("Cannot plot principal components having a 3- or higher dimensional domain.")
  
  # Wait for user input for each new eigenfunction
  oldPar <- graphics::par(no.readonly = TRUE)
  graphics::par(ask = TRUE)
  
  # set number of rows:
  # 1: all dimensions from left to right, "+" and "-" in one plot
  # 2: all dimensions from left to right, upper row for "+", lower for "-"
  nRows <- if(combined == TRUE)
  {
    if(all(dims == 1)) {1} else {
      warning("Cannot combine plots for two-dimensional elements. Will use separate plots (combined = FALSE)")
      2
    } 
  } else{2}
  
  graphics::par(mfrow = c(nRows, length(x$functions)))
  
  
  for(ord in plotPCs) # for each order
  {
    # calculate stretch factor if not given
    if(is.null(stretchFactor))
      stretchFactor <- stats::median(abs(x$scores[,ord]))
    
    PCplus <- x$meanFunction + stretchFactor * x$functions
    PCminus <- x$meanFunction - stretchFactor * x$functions
    
    for(rows in seq_len(nRows))
    {
      for(i in seq_len(length(x$functions))) # for each element
      {
        yRange <- range(PCplus[[i]]@X, PCminus[[i]]@X)
        main <- paste("PC", ord, "(explains", round(x$values[ord]/sum(x$values)*100, 2), "% of total variability)")
        
        if(dims[i] == 1)
        {
          funData::plot(x$meanFunction[[i]], lwd = 2, col = "black", 
                        main = main, ylim = yRange, ...)
          if(rows == 1)
            funData::plot(PCplus[[i]], obs = ord, 
                          add = TRUE, type = "p", pch = "+", col = "grey50", ...)
          if(rows == 2 | combined == TRUE)
            funData::plot( PCminus[[i]], obs = ord,
                           add = TRUE, type = "p", pch = "-", col = "grey50", ...)
        }  
        else # dims[i] == 2 (higher dimensional domains are not supported)
        {
          if(rows == 1)
            funData::plot(PCplus[[i]], obs = ord, main = main, ylim = yRange, ...)
          else
            funData::plot(PCminus[[i]], obs = ord, main = main, ylim = yRange, ...)
          
        }
      }
    }
  }
  
  graphics::par(oldPar)
  
  invisible(NULL)
}


#' Function prediction based on MFPCA results
#' 
#' Predict functions based on a truncated multivariate Karhunen-Loeve 
#' representation: \deqn{\hat x = \hat mu + \sum_{m = 1}^M \rho_m \hat \psi_m} 
#' with estimated mean function \eqn{\hat \mu} and principal components 
#' \eqn{\psi_m}. The scores \eqn{\rho_m} can be either estimated (reconstruction
#' of observed functions) or user-defined (construction of new functions).
#' 
#' @param object An object of class \code{MFPCAfit}, typically resulting from 
#'   a \link{MFPCA} function call.
#' @param scores A matrix containing the score values. The number of columns in 
#'   \code{scores} must equal the number of principal components in 
#'   \code{object}. Each row represents one curve. Defaults to the estimated
#'   scores in \code{object}, which yields reconstructions of the original
#'   data used for the MFPCA calculation.
#' @param ... Arguments passed to or from other methods.
#'   
#' @return A \code{multiFunData} object containing the predicted functions. 
#' 
#' @seealso \link{MFPCA}
#' 
#' @export
#' 
#' @examples 
#' #' # Simulate multivariate functional data on one-dimensonal domains
#' # and calculate MFPCA (cf. MFPCA help)
#' set.seed(1)
#' # simulate data (one-dimensional domains)
#' sim <-  simMultiFunData(type = "split", argvals = list(seq(0,1,0.01), seq(-0.5,0.5,0.02)),
#'                        M = 5, eFunType = "Poly", eValType = "linear", N = 100)
#' # MFPCA based on univariate FPCA
#' PCA <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"),
#'                                                      list(type = "uFPCA")))
#' 
#' # Reconstruct the original data
#' pred <- predict(PCA) # default reconstructs data used for the MFPCA fit
#' 
#' # plot the results: 1st element
#' plot(sim$simData[[1]]) # original data
#' plot(pred[[1]], add = TRUE, lty = 2) # reconstruction
#'
#' # plot the results: 2nd element
#' plot(sim$simData[[2]]) # original data
#' plot(pred[[2]], add = TRUE, lty = 2) # reconstruction
predict.MFPCAfit <- function(object, scores = object$scores, ...)
{
  if(!inherits(object, "MFPCAfit"))
    stop("Argument is not of class 'MFPCAfit'.")
  
  if(!(is.numeric(scores) & NCOL(scores) == length(object$values)))
     stop("Argument 'scores' must be a matrix with ", length(object$values), " columns.")
  
  return(object$meanFunction  + 
           multivExpansion(multiFuns = object$functions, scores = scores))
}


#' Print the results of a Multivariate Functional Principal Component Analysis
#' 
#' A \code{print} function for class \code{MFPCAfit}.
#' 
#' @param x An object of class \code{MFPCAfit}, usually returned by a
#'   call to \link{MFPCA}.
#' @param ... Arguments passed to or from other methods.
#'   
#' @export
print.MFPCAfit <- function(x, ...)
{
  if(!inherits(x, "MFPCAfit"))
    stop("Argument is not of class 'MFPCAfit'.")
  
  cat(funData::nObs(x$functions), "multivariate functional principal components estimated with",
      length(x$functions), "elements, each.\n", rep(c(" ", "*", " "), each = 10), "\n")
  
  cat("Eigenvalues:\n")
  print(x$values)
}


#' Summarize a Multivariate Functional Principal Component Analysis
#' 
#' A \code{summary} method for class \code{MFPCAfit}
#' 
#' @param object An object of class \code{MFPCAfit}, usually returned by a
#'   call to \link{MFPCA}.
#' @param ... Arguments passed to or from other methods.
#'   
#' @export 
#' @method summary MFPCAfit
summary.MFPCAfit <- function(object, ...)
{
  if(!inherits(object, "MFPCAfit"))
    stop("Argument is not of class 'MFPCAfit'.")
  
  vals <- object$values
  
  res <- rbind(`Eigenvalue` = vals, 
               `Proportion of variance explained` = vals/sum(vals),
               `Cumulative proportion` = cumsum(vals)/sum(vals))
  colnames(res) <- paste("PC", seq_len(length(vals)))
  
  attr(res, "npc") <- funData::nObs(object$functions)
  attr(res, "nel") <- length(object$functions)
  
  class(res) <- "summary.MFPCAfit"
  
  return(res)
}


#' Print summary of a Multivariate Functional Principal Component Analysis
#' 
#' A \code{print} method for class \code{MFPCAfit.summary}
#' 
#' @param x An object of class \code{MFPCAfit.summary}, usually returned by a
#'   call to \code{MFPCAfit.summary}.
#' @param ... Arguments passed to or from other methods.
#'   
#' @export
#' @method print summary.MFPCAfit
print.summary.MFPCAfit <- function(x, ...)
{
  if(!inherits(x, "summary.MFPCAfit"))
    stop("Argument is not of class 'summary.MFPCAfit'.")
  
  cat(attr(x, "npc"), "multivariate functional principal components estimated with",
      attr(x, "nel"), "elements, each.\n", rep(c(" ", "*", " "), each = 10), "\n")
  
  print.table(x, ...)
  
  invisible(x)
}