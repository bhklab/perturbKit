#' get_level5_ds - return the path to the level 5 dataset in a directory to simplify grabbing file names. 
#'
#' @param mypath Path to the file to which the downloaded data is saved
#'
#' @return path to level5 dataset in a directory. 
#' @export
#'

get_level5_ds <- function(mypath, mypattern=""){
  f <- list.files(mypath, pattern=sprintf("level5.*gctx"), ignore.case=TRUE)
  
  if (length(f) != 1){
    stop(sprintf("In get_level5_ds: Number of matched files in %s is %d, but should be 1", mypath, length(f)))
  } 
  return(file.path(mypath, f[1]))
}
