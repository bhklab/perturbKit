#' download_l1k_data - download data for L1000 from GEO and clue.io
#'
#' @param mypath Path to the file to which the downloaded data is saved
#' @param version numeric, either 2017 or 2020
#'
#' @export
#' @importFrom R.utils gunzip
download_l1k_data <- function(mypath, version=2017){
  if (getOption('timeout') < 100000){
    options(timeout=100000)
  }
  if (version == 2017){
    myfiles_exist <- (file.exists(file.path(mypath, "GSE92742_Broad_LINCS_sig_metrics.txt")) & 
                        file.exists(file.path(mypath, "GSE92742_Broad_LINCS_sig_info.txt")) &
                        file.exists(file.path(mypath, "GSE92742_Broad_LINCS_pert_metrics.txt")))
    
    if (!dir.exists(mypath)){
      dir.create(mypath)
    }
    if (!myfiles_exist){
      # Download L1000 files from the GEO repository GSE92742
      download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5FLevel5%5FCOMPZ%2EMODZ%5Fn473647x12328%2Egctx%2Egz", 
                    destfile=file.path(mypath, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz"))
      download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5FREADME%2Epdf", 
                    destfile=file.path(mypath, "GSE92742_Broad_LINCS_README.pdf"))
      download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fcell%5Finfo%2Etxt%2Egz", 
                    destfile=file.path(mypath, "GSE92742_Broad_LINCS_cell_info.txt.gz"))
      download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fgene%5Finfo%2Etxt%2Egz", 
                    destfile=file.path(mypath, "GSE92742_Broad_LINCS_gene_info.txt.gz"))
      download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fpert%5Finfo%2Etxt%2Egz", 
                    destfile=file.path(mypath, "GSE92742_Broad_LINCS_pert_info.txt.gz"))
      download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fpert%5Fmetrics%2Etxt%2Egz", 
                    destfile=file.path(mypath, "GSE92742_Broad_LINCS_pert_metrics.txt.gz"))
      download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fsig%5Finfo%2Etxt%2Egz", 
                    destfile=file.path(mypath, "GSE92742_Broad_LINCS_sig_info.txt.gz"))
      download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fsig%5Fmetrics%2Etxt%2Egz", 
                    destfile=file.path(mypath, "GSE92742_Broad_LINCS_sig_metrics.txt.gz"))
      download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5FLevel3%5FINF%5Fmlr12k%5Fn1319138x12328%2Egctx%2Egz", 
                    destfile=file.path(mypath, "GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz"))
      download.file(url="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5990023/bin/NIHMS917441-supplement-11.xlsx", 
                    destfile=file.path(mypath, "NIHMS917441-supplement-11.xlsx"))
      
      R.utils::gunzip(filename=paste(mypath, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz", sep="/"))
      R.utils::gunzip(filename=paste(mypath, "GSE92742_Broad_LINCS_cell_info.txt.gz", sep="/"))
      R.utils::gunzip(filename=paste(mypath, "GSE92742_Broad_LINCS_gene_info.txt.gz", sep="/"))
      R.utils::gunzip(filename=paste(mypath, "GSE92742_Broad_LINCS_pert_info.txt.gz", sep="/"))
      R.utils::gunzip(filename=paste(mypath, "GSE92742_Broad_LINCS_pert_metrics.txt.gz", sep="/"))
      R.utils::gunzip(filename=paste(mypath, "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="/"))
      R.utils::gunzip(filename=paste(mypath, "GSE92742_Broad_LINCS_sig_metrics.txt.gz", sep="/"))
      gunzip(filename=paste(mypath, "GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz", sep="/"))
    } else {
      print(sprintf("Data for version %s detected in %s, no download is necessary.", version, mypath))
    }
  } 
  else if (version == 2020){
    myfiles_exist <- (file.exists(file.path(mypath, "siginfo_beta.txt")) &
                        file.exists(file.path(mypath, "level5_beta_trt_cp_n720216x12328.gctx")))
    
    if (!dir.exists(mypath)){
      dir.create(mypath)
    }
    
    if (!myfiles_exist){
      # Paths extracted from clue.io
      download.file(url="https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/cellinfo_beta.txt", 
                    destfile=file.path(mypath, "cellinfo_beta.txt"))
      download.file(url="https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level5/level5_beta_trt_cp_n720216x12328.gctx", 
                    destfile=file.path(mypath, "level5_beta_trt_cp_n720216x12328.gctx"))
      download.file(url="https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/compoundinfo_beta.txt", 
                    destfile=file.path(mypath, "compoundinfo_beta.txt"))
      download.file(url="https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/geneinfo_beta.txt", 
                    destfile=file.path(mypath, "geneinfo_beta.txt"))
      download.file(url="https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/siginfo_beta.txt", 
                    destfile=file.path(mypath, "siginfo_beta.txt"))
      download.file(url="https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/instinfo_beta.txt", 
                    destfile=file.path(mypath, "instinfo_beta.txt"))
    } else {
      print(sprintf("Data for version %s detected in %s, no download is necessary.", version, mypath))
    }
  }
}


#' read_l1k_meta - load metadata for L1000 datasets
#'
#' @param mypath Path to data
#' @param version numeric, either 2017 or 2020
#'
#' @return list
#' @export
#'

read_l1k_meta <- function(mypath, version=2017){
  if (version == 2017){
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
                sigmetrics=sigmetrics, 
                landmarks=landmarks))
  } else if (version == 2020){
    cellinfo <- read.delim(file=file.path(mypath, "cellinfo_beta.txt"), stringsAsFactors = FALSE)
    geneinfo <- read.delim(file=file.path(mypath, "geneinfo_beta.txt"), stringsAsFactors = FALSE)
    pertinfo <- read.delim(file=file.path(mypath, "compoundinfo_beta.txt"), stringsAsFactors = FALSE)
    siginfo <- read.delim(file=file.path(mypath, "siginfo_beta.txt"), stringsAsFactors = FALSE)
    
    geneinfo$gene_id <- as.character(geneinfo$gene_id)
    # Remap geneinfo - because the 2020 data is beta, the column names may change.  I map them to the same ones as 2017 for consistency within Leo:
    colnames(geneinfo)[match(c("gene_id", "gene_symbol", "gene_title"), colnames(geneinfo))] <- c("pr_gene_id", "pr_gene_symbol", "pr_gene_title")
    geneinfo$pr_is_lm <- as.numeric(geneinfo$feature_space == "landmark")
    geneinfo$pr_is_bing <- as.numeric(geneinfo$feature_space %in% c("landmark", "best inferred"))
    
    # Remap siginfo, cellinfo, pertinfo
    colnames(siginfo)[match(c("cmap_name", "cell_iname"), colnames(siginfo))] <- c("pert_iname", "cell_id")
    colnames(cellinfo)[match(c("cell_iname"), colnames(cellinfo))] <- c("cell_id")
    colnames(pertinfo)[match(c("cmap_name"), colnames(cellinfo))] <- c("pert_iname")
    
    landmarks <- geneinfo[which(geneinfo$feature_space == "landmark"), ]
    
    return(list(cellinfo=cellinfo, 
                geneinfo=geneinfo, 
                pertinfo=pertinfo, 
                siginfo=siginfo,
                landmarks=landmarks))
  }
}
