
#' download the LINCS Kinomescan data
#' 
#' This function is disappointingly hacked, but I did not figure out a better way 
#' to extract the data from the LINCS HMS website. 
#' 
#' @param outdir Output directory
#' @export
#' 
#' @importFrom openxlsx read.xlsx

downloadKinomescan <- function(outdir="."){
  
  if (!dir.exists(outdir)){
    dir.create(outdir)
  }
  
  download.file("http://lincs.hms.harvard.edu/wordpress/wp-content/uploads/2013/11/HMS-LINCS_KinomeScan_Datasets_2018-01-18.xlsx",
                file.path(outdir, "HMS-LINCS_KinomeScan_Datasets_2018-01-18.xlsx"))
  
  
  kinomeManifest <- openxlsx::read.xlsx(file.path(outdir, "HMS-LINCS_KinomeScan_Datasets_2018-01-18.xlsx"))

  for (ii in seq_len(dim(kinomeManifest)[1])){
    
    download.file(sprintf("%s/results?search=&output_type=.csv", kinomeManifest$dataset_url[ii]), 
                  file.path(outdir, sprintf("kinomescan%s.csv", kinomeManifest$dataset_id[ii])))
    #https://lincs.hms.harvard.edu/db/datasets/20220/results?search=&output_type=.csv
  }
}


#' createKinomescanObj 
#' 
#' create the Kinomescan data object after downloading the data
#' 
#' @param datadir Data directory (the output directory from downloadKinomecsan)
#' @export





# This is a function to help track my data munging of the kinomescan data
scratchwork <- function(kinomedir=""){
  
  kinomef <- list.files(kinomedir, pattern = ".csv")
  
  kinomeData <- lapply(kinomef, FUN=function(x) read.csv(file.path(kinomedir, x)))
  
  # Dimension
  sapply(kinomeData, FUN=function(x) dim(x)[1])
  sapply(kinomeData, FUN=function(x) dim(x)[2])
}