#' read_l1k_meta
#' 
#' This function parses the L1000 metadata downloaded from either GEO 92742 or the clue.io data repository,
#' parses it, and returns the meta data as a list of dataframes. 
#' @param mypath A string with the path to the directory containing the metadata, with subfolders 2017 and 2020.
#' @param version A numeric indicating which version of the data to use - 2017 or 2020 (defaults to 2017).

read_l1k_meta <- function(mypath, version=2017){
  if (version == 2017){
    #mypath <- file.path(mypath, "2017")
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
    #mypath <- file.path(mypath, "2020")
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