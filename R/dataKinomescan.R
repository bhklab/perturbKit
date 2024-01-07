
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

createKinomescanObj <- function(datadir, outdir="."){
  
  kinomef <- list.files(kinomedir, pattern = ".csv")
  
  kinomeData <- lapply(kinomef, FUN=function(x) read.csv(file.path(kinomedir, x)))
  
  # The kinomescan data is organized with three distinct sets of compound assays:
  # One set with 7 columns that measures %binding compared to control
  # One set with 6 columns that measures Kd. 
  #
  ix <- which(sapply(kinomeData, FUN=function(x) "X..Control" %in% colnames(x)))
  jx <- which(sapply(kinomeData, FUN=function(x) "Kd" %in% colnames(x)))
  
  # Get all proteins against which the kinomescan drugs were assayed
  myprots <- purrr::reduce(lapply(kinomeData[ix], FUN=function(x) x[,c('Protein.HMS.LINCS.ID', 'Protein.Name')]), rbind)
  
  protset <- dplyr::distinct(myprots)
  drugData <- data.frame(SmallMoleculeHMSLINCSID=character(), 
                         SmallMoleculeName=character(),
                         uMDose=character())
  pctControlData <- numeric()
  
  for (ii in seq_along(kinomeData)){
    # Only consider kinomescan results with 7 columns, i.e. containing X..Control
    if (dim(kinomeData[[ii]])[2] == 7){
      
      print(ii)
      mydoses <- unique(kinomeData[[ii]]$Assay.compound.conc)
      
      for (jj in mydoses){
        kx <- which(kinomeData[[ii]]$Assay.compound.conc == jj)
        
        sub <- kinomeData[[ii]][kx,]
        
        if (dim(dplyr::distinct(sub[, c("Small.Molecule.HMS.LINCS.ID", "Small.Molecule.Name", 
                                        "Assay.compound.conc", "Conc.unit")]))[1] != 1){
          stop(sprintf("Kinome Data subset does not have unique compound values, ii = %d", ii))
        }
        
        subDoseuM <- sub$Assay.compound.conc[1]
        if (sub$Conc.unit[1] == "nM"){
          subDoseuM <- subDoseuM * 1000
        }
        
        drugData <- rbind(drugData, data.frame(SmallMoleculeHMSLINCSID = sub$Small.Molecule.HMS.LINCS.ID[1],
                                               SmallMoleculeName = sub$Small.Molecule.Name[1],
                                               uMDose = subDoseuM))
        pctControlData <- cbind(pctControlData, 
                                sub[match(protset$Protein.HMS.LINCS.ID, sub$Protein.HMS.LINCS.ID), "X..Control"])
        
      }
    }
  }
  
  
  # Get dataset with kD measurements
  kdp <- purrr::reduce(lapply(kinomeData[jx], FUN=function(x) x[,c('Protein.HMS.LINCS.ID', 'Protein.Name')]), rbind)
  kdprotset <- dplyr::distinct(kdp)
  
  kdDrugData <- data.frame(SmallMoleculeHMSLINCSID=character(), 
                           SmallMoleculeName=character())
  kdData <- numeric()
  
  for (ii in seq_along(kinomeData)){
    # Only consider kinomescan results with 6 columns, i.e. containing Kd
    if (dim(kinomeData[[ii]])[2] == 6){
      
      print(ii)
      sub <- kinomeData[[ii]]
        
      if (dim(dplyr::distinct(sub[, c("Small.Molecule.HMS.LINCS.ID", "Small.Molecule.Name")]))[1] != 1){
        stop(sprintf("Kinome Data subset does not have unique compound values, ii = %d", ii))
      }
        
      kdDrugData <- rbind(kdDrugData, data.frame(SmallMoleculeHMSLINCSID = sub$Small.Molecule.HMS.LINCS.ID[1],
                                               SmallMoleculeName = sub$Small.Molecule.Name[1]))
      kdData <- cbind(kdData, sub[match(kdprotset$Protein.HMS.LINCS.ID, sub$Protein.HMS.LINCS.ID), "Kd"])
    
    }
  }
  
}




# This is a function to help track my data munging of the kinomescan data
scratchwork <- function(kinomedir=""){
  
  kinomef <- list.files(kinomedir, pattern = ".csv")
  
  kinomeData <- lapply(kinomef, FUN=function(x) read.csv(file.path(kinomedir, x)))
  
  # Dimension
  sapply(kinomeData, FUN=function(x) dim(x)[1])
  sapply(kinomeData, FUN=function(x) dim(x)[2])
}