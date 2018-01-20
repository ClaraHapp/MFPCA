#' Plot the Scores of a Multivariate Functional Principal Component Analysis
#' 
#' This function plots two scores of a multivariate functional principal 
#' component analysis for each observation.
#' 
#' @param MFPCAobj An object of class \code{MFPCAfit}, typically returned by the
#'   \link{MFPCA} function.
#' @param choices The indices of the scores that should by displayed. Defaults 
#'   to \code{1:2}, i.e. the scores corresponding to the two leading modes of 
#'   variability in the data.
#' @param scale Logical. Should the scores be scaled by the estimated
#'   eigenvalues to emphasize the proportions of total variance explained by the
#'   components. Defaults to \code{FALSE}.
#'   
#' @return A bivariate plot of scores.
#' 
#' @seealso \code{\link{MFPCA}}
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
#' scoreplot(PCA, scaling = TRUE) # scale the scores by the first two eigenvalues 
scoreplot.MFPCAfit <- function(MFPCAobj, choices = 1:2, scale = FALSE, ...)
{
  if(length(choices) != 2)
    stop("Length of choices must be 2.")
  
  if(NCOL(MFPCAobj$scores) < max(choices))
    stop(paste("Argument choices requires", max(choices), 
               "scores, MFPCA object contains only", NCOL(MFPCAobj$scores), "."))
  
  # check for labels, otherwise construct them TODO: rownames for scores in MFPCA
  lab <- {if(is.null(rownames(MFPCAobj$scores)))
    1:NROW(MFPCAobj$scores)
    else
      rownames(MFPCAobj$scores)}
  
  # scale scores if scale = TRUE
  plotScore <- {
    if(scale) # multiply row-wise by eigenvalues
      t(MFPCAobj$values[choices] * t(MFPCAobj$scores[, choices]))
    else
      MFPCAobj$scores[, choices]
  }
  
  # plot scores
  graphics::plot.default(x = plotScore, y = NULL, type = "n", 
                         xlab = paste("Scores PC", choices[1]), ylab = paste("Scores PC", choices[2]), ...)
  graphics::text(x = plotScore, y = NULL, labels = lab, ...)
  
  abline(h = 0, lty = 3, col = "gray")
  abline(v = 0, lty = 3, col = "gray")
  
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
#' @param MFPCAobj An object of class \code{MFPCAfit}, typically returned by the
#'   \link{MFPCA} function.
#' @param plotPCs The principal components to be plotted. Defaults to all 
#'   components in the \code{MFPCAfit} object.
#' @param stretchFactor The factor by which the principal components are 
#'   multiplied before adding / subtracting them from the mean function. If 
#'   \code{NULL} (the default), the median absolute value of the scores of each 
#'   eigenfunction is used.
#' @param combined Logical: Should the plots be combined? (Works only if all 
#'   dimensions are one-dimensional). Defaults to \code{FALSE}.
#' @param .. Further graphical parameters passed to the \link[funData]{plot}
#'   functions for functional data.
#'   
#' @return A plot of the principal components as perturbations of the mean.
#' 
#' @seealso \link{MFPCA}, \link[funData]{plot}
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
plot.MFPCAfit <- function(MFPCAobj, plotPCs = 1:nObs(MFPCAobj$functions), stretchFactor = NULL, combined = FALSE, ...)
{
  # check dimensions
  dims <- dimSupp(MFPCAobj$functions)
  
  if(any(dims > 2))
    stop("Cannot plot principal components having a 3- or higher dimensional domain.")
  
  # Wait for user input for each new eigenfunction
  oldPar <- par(no.readonly = TRUE)
  par(ask = TRUE)
  
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
  
  par(mfrow = c(nRows, length(MFPCAobj$functions)))
  
  
  for(ord in plotPCs) # for each order
  {
    # calculate stretch factor if not given
    if(is.null(stretchFactor))
      stretchFactor <- median(abs(MFPCAobj$scores[,ord]))
    
    PCplus <- MFPCAobj$meanFunction + stretchFactor * MFPCAobj$functions
    PCminus <- MFPCAobj$meanFunction - stretchFactor * MFPCAobj$functions
    
    for(rows in 1:nRows)
    {
      for(i in 1:length(MFPCAobj$functions)) # for each element
      {
        yRange <- range(PCplus[[i]]@X, PCminus[[i]]@X)
        main <- paste("PC", ord, "(explains", round(MFPCAobj$values[ord]/sum(MFPCAobj$values)*100, 2), "% of total variability)")
        
        if(dims[i] == 1)
        {
          plot(MFPCAobj$meanFunction[[i]], lwd = 2, col = "black", 
               main = main, ylim = yRange, ...)
          if(rows == 1)
            plot(PCplus[[i]], obs = ord, 
                 add = TRUE, type = "p", pch = "+", col = "grey50", ...)
          if(rows == 2 | combined == TRUE)
            plot( PCminus[[i]], obs = ord,
                  add = TRUE, type = "p", pch = "-", col = "grey50", ...)
        }  
        else # dims[i] == 2 (higher dimensional domains are not supported)
        {
          if(rows == 1)
            plot(PCplus[[i]], obs = ord, main = main, ylim = yRange, ...)
          else
            plot(PCminus[[i]], obs = ord, main = main, ylim = yRange, ...)
          
        }
      }
    }
  }
  
  par(oldPar)
  
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
#' @param MFPCAobj An object of class \code{MFPCAfit}, typically resulting from 
#'   a \link{MFPCA} function call.
#' @param scores A matrix containing the score values. The number of columns in 
#'   \code{scores} must equal the number of principal components in 
#'   \code{MFPCAobj}. Each row represents one curve. Defaults to the estimated
#'   scores in \code{MFPCAobj}, which yields reconstructions of the original
#'   data used for the MFPCA calculation.
#'   
#' @return A \code{multiFunData} object containing the predicted functions. 
#' 
#' @seealso \link{MFPCA}
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
predict.MFPCAfit <- function(MFPCAobj, scores = MFPCAobj$scores)
{
  return(MFPCAobj$meanFunction  + 
           multivExpansion(multiFuns = MFPCAobj$functions, scores = scores))
}


#' Print the results of a Multivariate Functional Principal Component Analysis
#' 
#' A \code{print} function for class \code{MFPCAfit}.
#' 
#' @param MFPCAobj An object of class \code{MFPCAfit}, usually returned by a
#'   call to \link{MFPCA}.
print.MFPCAfit <- function(MFPCAobj)
{
  cat(nObs(MFPCAobj$functions), "multivariate functional principal components estimated with",
      length(MFPCAobj$functions), "elements, each.\n", rep(c(" ", "*", " "), each = 10), "\n")
  
  cat("Eigenvalues:\n")
  print(MFPCAobj$values)
}


#' Summarize a Multivariate Functional Principal Component Analysis
#' 
#' A \code{summary} method for class \code{MFPCAfit}
#' 
#' @param MFPCAobj An object of class \code{MFPCAfit}, usually returned by a
#'   call to \link{MFPCA}.
summary.MFPCAfit <- function(MFPCAobj)
{
  cat(nObs(MFPCAobj$functions), "multivariate functional principal components estimated with",
      length(MFPCAobj$functions), "elements, each.\n", rep(c(" ", "*", " "), each = 10), "\n")
  
  vals <- MFPCAobj$values
  
  res <- rbind(`Eigenvalue` = vals, 
               `Proportion of variance explained` = vals/sum(vals),
               `Cumulative proportion` = cumsum(vals)/sum(vals))
  colnames(res) <- paste("PC", 1:length(vals))
  
  return(res)
}