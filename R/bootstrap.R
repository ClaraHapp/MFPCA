#### Functions for bootstrapping in MFPCA ####

#' Sample stratified indices according to a factor variable
#' 
#' @param f A factor variable
#'   
#' @return A vector of the same length as \code{f}, containing the resampled
#'   indices, stratified according to the levels of f
#'   
#' @keywords internal
#'   
#' @examples
#' 
#' # create factor variable
#' f <- as.factor(c(1,1,1,2,4,4,4,4,6,6))
#' table(f)
#' 
#' sampleInd <- MFPCA:::stratSample(f)
#' table(f[sampleInd])
stratSample <- function(f)
{
  if(!is.factor(f))
    stop("f has to be a factor variable.")
  
  return( unlist(lapply(split(1:length(f), f), 
                        function(i){i[sample.int(length(i), replace = TRUE)]}), # use this in order to avoid surprises if length(i) = 1 (see sample help)
                 use.names = FALSE) )
}