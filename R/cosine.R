#' cosine
#' 
#' This function computes the cosine distance between all pairs of columns of two input matrices
#' @param mat1 Numeric, a matrix of size N x d1
#' @param mat2 Numeric, a matrix of size N x d2
#'
#' @return matrix
#'
#' @export
#' 

cosine <- function(mat1, mat2){
  a <- t(mat1) %*% mat2
  
  mag1 <- sqrt(diag(t(mat1) %*% mat1))
  mag2 <- sqrt(diag(t(mat2) %*% mat2))
  
  return(a/outer(mag1, mag2, "*"))
}