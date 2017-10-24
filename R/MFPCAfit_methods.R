biplot.MFPCAfit <- function(MFPCAobj, choices = 1:2, ...)
{
  if(length(choices) != 2)
    stop("Length of choices must be 2.")
  
  if(NCOL(MFPCAobj$scores) < max(choices))
    stop(paste("Argument choices requires", max(choices), 
               "scores, MFPCA object contains only", NCOL(MFPCAobj$scores), "."))
  
  graphics:::plot.default(x = MFPCAobj$scores[, choices], y = NULL, ...)
    
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

