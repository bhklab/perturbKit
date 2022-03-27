#' get_top_k
#' 
#' This is a helper function used in compute_sim_block that identifies the row identifiers of the K largest or smallest elements in each column of a matrix.
#'
#' @param dsmat A matrix
#' @param k A non-negative integer indicating the number of largest elements to return
#' @param decreasing A boolean indicating whether to get the largest elements (TRUE) or the smallest elements (FALSE); true by default. 
#' 
#' @return vector
#' @export

get_top_k <- function(dsmat, k=50, decreasing=TRUE){
  
  if (k > 500){
    k <- 500
  }
  mysets <- lapply(seq(dim(dsmat)[2]), FUN=function(x) rownames(dsmat)[order(dsmat[,x], decreasing=decreasing)[seq(k)]])
  return(mysets)
}
