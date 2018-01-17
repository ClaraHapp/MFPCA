biplot.MFPCAfit <- function(MFPCAobj, choices = 1:2, ...)
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
  
  
  graphics::plot.default(x = MFPCAobj$scores[, choices], y = NULL, type = "n", 
                         xlab = paste("Scores PC", choices[1]), ylab = paste("Scores PC", choices[2]), ...)
  graphics::text(x = MFPCAobj$scores[, choices], y = NULL, labels = lab, ...)
  
  abline(h = 0, lty = 3, col = "gray")
  abline(v = 0, lty = 3, col = "gray")
  
  return(invisible)
}

#' @param plotPCs The principal components to be plotted. By default all components in the \code{MFPCAfit} object
#' @param stretchFactor The factor by which the principal components are 
#'   multiplied before adding / subtracting them from the mean function. If
#'   \code{NULL} (the default), the median absolute value of the scores is used.
plot.MFPCAfit <- function(MFPCAobj, plotPCs = 1:nObs(MFPCAobj$functions), stretchFactor = NULL , ...)
{
  dims <- dimSupp(MFPCAobj$functions)
  if(any(dims > 2))
    stop("Cannot plot principal components having a 3- or higher dimensional domain.")
  else
    allDims1 <- all(dims == 1)
  oldPar <- par(no.readonly = TRUE)
  
  par(ask = TRUE)
  
  if(allDims1)
  {
    # all dimensions from left to right
    par(mfrow = c(1, length(MFPCAobj$functions)))
  }
  else
  {
    # all dimensions from left to right, upper row for "+", lower for "-"
    par(mfrow = c(2, length(MFPCAobj$functions)))
  }
  
  for(ord in plotPCs) # for each order
  {
    # calculate stretch factor if not given
    if(is.null(stretchFactor))
      stretchFactor <- median(abs(MFPCAobj$scores[,ord]))
    
    PCplus <- MFPCAobj$meanFunction + stretchFactor * MFPCAobj$functions
    PCminus <- MFPCAobj$meanFunction - stretchFactor * MFPCAobj$functions
    
    for(rows in 1:ifelse(allDims1 == TRUE, 1, 2))
    {
      for(i in 1:length(MFPCAobj$functions)) # for each element
      {
        if(dims[i] == 1)
        {
          plot(MFPCAobj$meanFunction[[i]], lwd = 2, col = "black", 
               main = paste("PC", ord, "(explains", round(MFPCAobj$values[ord]/sum(MFPCAobj$values)*100, 2), "% of total variability)"),
               ylim =  range(PCplus[[i]]@X, PCminus[[i]]@X), ...)
          if(rows == 1)
            plot(PCplus[[i]], obs = ord, 
                 add = TRUE, type = "p", pch = "+", col = "grey50", ...)
          if(allDims1 == TRUE | rows == 2)
            plot( PCminus[[i]], obs = ord,
                  add = TRUE, type = "p", pch = "-", col = "grey50", ...)
        }  
        else # dims[i] == 2 (higher dimensional domains are not supported)
        {
          if(rows == 1)
            plot(PCplus[[i]], obs = ord, main = paste("PC", ord, "(explains", round(MFPCAobj$values[ord]/sum(MFPCAobj$values)*100, 2), "% of total variability)"),
                 ylim = range(PCplus[[i]]@X, PCminus[[i]]@X), ...)
          else
            plot(PCminus[[i]], obs = ord, main = paste("PC", ord, "(explains", round(MFPCAobj$values[ord]/sum(MFPCAobj$values)*100, 2), "% of total variability)"),
                 ylim = range(PCplus[[i]]@X, PCminus[[i]]@X), ...)
          
        }
      }
    }
  }
  
  par(oldPar)
}

predict.MFPCAfit <- function(MFPCAobj, scores = MFPCAobj$scores)
{
  return(MFPCAobj$meanFunction  + 
           multivExpansion(multiFuns = MFPCAobj$functions, scores = scores))
}

print.MFPCAfit <- function(MFPCAobj)
{
  cat(nObs(MFPCAobj$functions), "multivariate functional principal components estimated with",
      length(MFPCAobj$functions), "elements, each.\n", rep(c(" ", "*", " "), each = 10), "\n")
  
  cat("Eigenvalues:\n")
  print(MFPCAobj$values)
}

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

