#' Calculate multivariate basis expansion
#' 
#' @param multiFuns A multivariate functional data object, containing the
#'   multivariate basis functions
#' @param scores A matrix containing the scores for each observation in each
#'   row. The number of columns must match the number of basis functions.
#'   
#' @return A \code{multiFunData} object containing the expanded functions for
#'   each observation.
multivExpansion <- function(multiFuns, scores)
{
  if(nObs(multiFuns) != NCOL(scores))
    stop("Number of scores does not match number of eigenfunctions.")
  
  # calculate linear combination of multivariate basis functions
  univExp <- foreach::foreach(j = 1:length(multiFuns)) %do% { # %do% might require extra loading
    univExpansion(type = "default", 
                  scores = scores,
                  functions = multiFuns[[j]])
  }
  
  # return as multiFunData object
  return(multiFunData(univExp))
}