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

plot.MFPCAfit <- function()
{
  
  
  
}

predict.MFPCAfit <- function()
{
  
  
  
}

print.MFPCAfit <- function()
{
  
  
  
}

summary.MFPCAfit <- function()
{
  
  
  
}

