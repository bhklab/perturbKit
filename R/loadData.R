download_l1k_data <- function(mypath){
  myfiles_exist <- (file.exists(paste(mypath, "GSE92742_Broad_LINCS_sig_metrics.txt", sep="/")) & 
                      file.exists(paste(mypath, "GSE92742_Broad_LINCS_sig_info.txt", sep="/")) &
                      file.exists(paste(mypath, "GSE92742_Broad_LINCS_pert_metrics.txt", sep="/")))
  if (!myfiles_exist){
    # Download L1000 files from the GEO repository GSE92742
    download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5FLevel5%5FCOMPZ%2EMODZ%5Fn473647x12328%2Egctx%2Egz", 
                  destfile=paste(mypath, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz", sep="/"))
    download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5FLevel3%5FINF%5Fmlr12k%5Fn1319138x12328%2Egctx%2Egz", 
                  destfile=paste(mypath, "GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz", sep="/"))
    
    download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5FREADME%2Epdf", 
                  destfile=paste(mypath, "GSE92742_Broad_LINCS_README.pdf", sep="/"))
    download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fcell%5Finfo%2Etxt%2Egz", 
                  destfile=paste(mypath, "GSE92742_Broad_LINCS_cell_info.txt.gz", sep="/"))
    download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fgene%5Finfo%2Etxt%2Egz", 
                  destfile=paste(mypath, "GSE92742_Broad_LINCS_gene_info.txt.gz", sep="/"))
    download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fpert%5Finfo%2Etxt%2Egz", 
                  destfile=paste(mypath, "GSE92742_Broad_LINCS_pert_info.txt.gz", sep="/"))
    download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fpert%5Fmetrics%2Etxt%2Egz", 
                  destfile=paste(mypath, "GSE92742_Broad_LINCS_pert_metrics.txt.gz", sep="/"))
    download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fsig%5Finfo%2Etxt%2Egz", 
                  destfile=paste(mypath, "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="/"))
    download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fsig%5Fmetrics%2Etxt%2Egz", 
                  destfile=paste(mypath, "GSE92742_Broad_LINCS_sig_metrics.txt.gz", sep="/"))
    download.file(url="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5990023/bin/NIHMS917441-supplement-11.xlsx", 
                  destfile=paste(mypath, "NIHMS917441-supplement-11.xlsx", sep="/"))
    
    gunzip(destfile=paste(mypath, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz", sep="/"))
    gunzip(destfile=paste(mypath, "GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz", sep="/"))
    
    gunzip(destfile=paste(mypath, "GSE92742_Broad_LINCS_cell_info.txt.gz", sep="/"))
    gunzip(destfile=paste(mypath, "GSE92742_Broad_LINCS_gene_info.txt.gz", sep="/"))
    gunzip(destfile=paste(mypath, "GSE92742_Broad_LINCS_pert_info.txt.gz", sep="/"))
    gunzip(destfile=paste(mypath, "GSE92742_Broad_LINCS_pert_metrics.txt.gz", sep="/"))
    gunzip(destfile=paste(mypath, "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="/"))
    gunzip(destfile=paste(mypath, "GSE92742_Broad_LINCS_sig_metrics.txt.gz", sep="/"))
    
  }
}


read_l1k_meta <- function(mypath){
  cellinfo <- read.delim(file=paste(mypath, "GSE92742_Broad_LINCS_cell_info.txt", sep="/"), stringsAsFactors=FALSE)
  geneinfo <- read.delim(file=paste(mypath, "GSE92742_Broad_LINCS_gene_info.txt", sep="/"), stringsAsFactors=FALSE)
  pertinfo <- read.delim(file=paste(mypath, "GSE92742_Broad_LINCS_pert_info.txt", sep="/"), stringsAsFactors=FALSE)
  siginfo <- read.delim(file=paste(mypath, "GSE92742_Broad_LINCS_sig_info.txt", sep="/"), stringsAsFactors=FALSE)
  instinfo <- read.delim(file=paste(mypath, "GSE92742_Broad_LINCS_inst_info.txt", sep="/"), stringsAsFactors=FALSE)
  
  geneinfo$pr_gene_id <- as.character(geneinfo$pr_gene_id)
  
  pertmetrics <- read.delim(file=paste(mypath, "GSE92742_Broad_LINCS_pert_metrics.txt", sep="/"), stringsAsFactors=FALSE)
  sigmetrics <- read.delim(file=paste(mypath, "GSE92742_Broad_LINCS_sig_metrics.txt", sep="/"), stringsAsFactors=FALSE)
  
  landmarks <- geneinfo[which(geneinfo$pr_is_lm == 1), ]
  
  return(list(cellinfo=cellinfo, 
              geneinfo=geneinfo,
              pertinfo=pertinfo,
              siginfo=siginfo,
              instinfo=instinfo,
              pertmetrics=pertmetrics,
              sigmetrics=sigmetrics))
}