#' calcSimBlock
#' 
#' This function takes as input two gct datasets and computes a similarity matrix between
#' columns of the first dataset and columns of the second.  It requires the gene space to be 
#' common to both.  Note that some similarity metrics, like GSEA-based metrics, are not symmetric.
#' ds1 and ds2 may also be matrices of dimension N_features x N_signatures. 
#' @param ds1 A gct instance, dataset 1 or a matrix
#' @param ds2 A gct instance, dataset 2 (could be the same as dataset 1) or a matrix (must be same datatype as ds1)
#' @param metric A string, one of wtcs, cosine, spearman, pearson, fastwtcs
#' @param kgenes For those metrics that restrict the number of genes (wtcs, cosine), a parameter that determines how many genes are used to compute similarity. Default is 50 for GSEA, 0 (all) for cosine.
#' @param gseaParam The GSEA parameter that affects the weighting of the gene values in GSEA.  Default is 1.
#' @param nperms Numeric, indicates how many permutations are run by fgsea for significance calculations.  As p-values from GSEA are not used, this defaults to 1.
#' @param parallel Logical, whether to parallelize. Default = 0.
#' @param numCores numeric. If parallel, how many cores to use. Default = parallel::detectCores()
#'
#' @return matrix
#'
#' @export
#' @importFrom coop cosine
#' @importFrom parallel mclapply detectCores

calcSimBlock <- function(ds1, ds2, metric="cosine", kgenes=0, gseaParam=1, nperms=1, parallel=0, numCores=parallel::detectCores()){

  metric <- match.arg(metric, c("wtcs", "cosine", "spearman", "pearson", "fastwtcs"))
  
  # Consider deprecating parallelization; fgsea doesn't seem to benefit from it, e.g. It's far
  # better for the user to parallelize any matrices they want to compute rather than trying 
  # to compute a single block. Perhaps a wrapper could parallelize matrix computation by dividing
  # the problem into subproblems. 
  
  # This is a bit of a hack, but it casts input matrices into gct objects. Probably better to rewrite
  # the function to operate on matrices, and extract matrices from GCT objects. 
  if (is.numeric(ds1) & is.numeric(ds2)){
    ds1 <- cmapR::GCT(mat=ds1, rid=as.character(seq(dim(ds1)[1])), cid=as.character(seq(dim(ds1)[2])))
    ds2 <- cmapR::GCT(mat=ds2, rid=as.character(seq(dim(ds2)[1])), cid=as.character(seq(dim(ds2)[2])))
  }
  
  if (!all.equal(ds1@rid, ds2@rid)){
    return("Error: input datasets did not have the same row space or rids.")
  }
  
  # Replace row names with integers for memory purposes
  rownames(ds1@mat) <- seq_len(length(ds1@rid))
  rownames(ds2@mat) <- seq_len(length(ds2@rid))
    
  if (parallel == 0){
    if (metric == "slowcosine"){
      
      if (kgenes > 0){
        ds1@mat <- as.matrix(ds1@mat * (MatrixGenerics::colRanks(ds1@mat, preserveShape = TRUE) <= kgenes | 
                                        MatrixGenerics::colRanks(ds1@mat, preserveShape = TRUE) > (dim(ds1@mat)[1] - kgenes)))
      }
      
      cosres <- matrix(numeric(dim(ds1@mat)[2] * dim(ds2@mat)[2]), nrow = dim(ds1@mat)[2], dimnames=list(ds1@cid, ds2@cid))
      # Currently iterative; parallelize
      for (ii in seq(dim(ds1@mat)[2])){
        cosres[ii,] <- sapply(seq(dim(ds2@mat)[2]), FUN=function(x) coop::cosine(ds1@mat[,ii], ds2@mat[,x]))
      }
      cosres <- t(cosres)
      colnames(cosres) <- ds1@cid
      rownames(cosres) <- ds2@cid
      return(cosres)
      
    } else if (metric == "cosine"){
      
      if (kgenes > 0){
        ds1@mat <- as.matrix(ds1@mat * (MatrixGenerics::colRanks(ds1@mat, preserveShape = TRUE) <= kgenes | 
                                        MatrixGenerics::colRanks(ds1@mat, preserveShape = TRUE) > (dim(ds1@mat)[1] - kgenes)))
      }
      
      cosres <- perturbKit::cosine(ds1@mat, ds2@mat)
      colnames(cosres) <- ds1@cid
      rownames(cosres) <- ds2@cid
      return(cosres)

    } else if (metric == "wtcs"){
      
      upsets <- getTopK(ds1@mat, k=kgenes, decreasing=TRUE)
      dnsets <- getTopK(ds1@mat, k=kgenes, decreasing=FALSE)
      allsets <- c(upsets, dnsets)
      names(allsets) <- as.character(seq(length(allsets)))
      ncol <- dim(ds1@mat)[2]
      
      es <- matrix(numeric(dim(ds2@mat)[2] * ncol), ncol = ncol)
      
      for (ii in seq(dim(ds2@mat)[2])){
        tempGSEA <- fgsea(allsets, ds2@mat[, ii], nperm=nperms, gseaParam=gseaParam)
        # This generates p-values, but is more expensive. 
        # tGSEA <- fgseaMultilevel(allsets, ds2@mat[, ii], sampleSize=k)
        upscore <- tempGSEA$ES[1:ncol]
        dnscore <- tempGSEA$ES[(1+ncol):(2*ncol)]
        es[ii, ] <- (upscore - dnscore)/2 * abs(sign(upscore) - sign(dnscore))/2
        
        if (ii %% 100 == 0){
          print(sprintf("%d / %d", ii, dim(ds2@mat)[2]))
        }
      }
      
      rownames(es) <- ds2@cid
      colnames(es) <- ds1@cid
      return(es)
      
    } else if (metric == "fastwtcs"){
      upsets <- getTopK(ds1@mat, k=kgenes, decreasing=TRUE)
      dnsets <- getTopK(ds1@mat, k=kgenes, decreasing=FALSE)
      allsets <- c(upsets, dnsets)
      names(allsets) <- as.character(seq(length(allsets)))
      ncol <- dim(ds1@mat)[2]
      
      es <- matrix(numeric(dim(ds2@mat)[2] * ncol), ncol = ncol)
      #print(sprintf("ds1 mat dimension: %d x %d \n", dim(ds1@mat)[1], dim(ds1@mat)[2]))
      #print(sprintf("ds2 mat dimension: %d x %d \n", dim(ds2@mat)[1], dim(ds2@mat)[2]))
      
      for (ii in seq(dim(ds2@mat)[2])){
        tempGSEA <- fgseaBase(allsets, ds2@mat[, ii], nperm=nperms, gseaParam=gseaParam)
        # This generates p-values, but is more expensive. 
        # tGSEA <- fgseaMultilevel(allsets, ds2@mat[, ii], sampleSize=k)
        upscore <- tempGSEA$ES[1:ncol]
        dnscore <- tempGSEA$ES[(1+ncol):(2*ncol)]
        es[ii, ] <- (upscore - dnscore)/2 * abs(sign(upscore) - sign(dnscore))/2
        
        if (ii %% 100 == 0){
          print(sprintf("%d / %d", ii, dim(ds2@mat)[2]))
        }
      }
      
      rownames(es) <- ds2@cid
      colnames(es) <- ds1@cid
      return(es)
    } else if (metric == "pearson") {
      
      pearson_res <- cor(ds1@mat, ds2@mat, method="pearson")
      pearson_res <- t(pearson_res)
      colnames(pearson_res) <- ds1@cid
      rownames(pearson_res) <- ds2@cid
      return(pearson_res)
      
    } else if (metric == "spearman") {
      
      spearman_res <- cor(ds1@mat, ds2@mat, method="spearman")
      spearman_res <- t(spearman_res)
      colnames(spearman_res) <- ds1@cid
      rownames(spearman_res) <- ds2@cid
      return(spearman_res)
      
    }
  } 
  else {
    # Parallelized
    
    if (metric == "slowcosine"){
      
      if (kgenes > 0){
        ds1@mat <- as.matrix(ds1@mat * (MatrixGenerics::colRanks(ds1@mat, preserveShape = TRUE) <= kgenes | 
                                        MatrixGenerics::colRanks(ds1@mat, preserveShape = TRUE) > (dim(ds1@mat)[1] - kgenes)))
      }
      
      x <- parallel::mclapply(seq(dim(ds1@mat)[2]), FUN=function(x) {
        sapply(seq(dim(ds2@mat)[2]), FUN=function(y) coop::cosine(ds1@mat[,x], ds2@mat[,y]))
        }, mc.cores=numCores)
      cosres <- matrix(unlist(x), nrow=dim(ds1@mat)[2])
      
      cosres <- t(cosres)
      colnames(cosres) <- ds1@cid
      rownames(cosres) <- ds2@cid
      return(cosres)
      
    } else if (metric == "cosine"){
      #lol no 
      
      return(calcSimBlock(ds1, ds2, metric="cosine", kgenes=kgenes, parallel=0))
      
    } else if (metric == "wtcs"){
      
      upsets <- getTopK(ds1@mat, k=kgenes, decreasing=TRUE)
      dnsets <- getTopK(ds1@mat, k=kgenes, decreasing=FALSE)
      allsets <- c(upsets, dnsets)
      names(allsets) <- as.character(seq(length(allsets)))
      ncol <- dim(ds1@mat)[2]

      tempGSEA <- parallel::mclapply(seq(dim(ds2@mat)[2]), FUN=function(x) fgsea(allsets, ds2@mat[,x], nperm=nperms, gseaParam=gseaParam), mc.cores=numCores)
      zes <- lapply(tempGSEA, FUN=function(y) (y$ES[1:ncol] - y$ES[(1+ncol):(2*ncol)])/2 * abs(sign(y$ES[1:ncol]) - sign(y$ES[(1+ncol):(2*ncol)]))/2)
      es <- t(matrix(unlist(zes), ncol = dim(ds2@mat)[2]))
      
      rownames(es) <- ds2@cid
      colnames(es) <- ds1@cid
      return(es)
      
    } else if (metric == "fastwtcs"){
      upsets <- getTopK(ds1@mat, k=kgenes, decreasing=TRUE)
      dnsets <- getTopK(ds1@mat, k=kgenes, decreasing=FALSE)
      allsets <- c(upsets, dnsets)
      names(allsets) <- as.character(seq(length(allsets)))
      ncol <- dim(ds1@mat)[2]
      
      tempGSEA <- parallel::mclapply(seq(dim(ds2@mat)[2]), FUN=function(x) fgseaBase(allsets, ds2@mat[,x], nperm=nperms, gseaParam=gseaParam), mc.cores=numCores)
      zes <- lapply(tempGSEA, FUN=function(y) (y$ES[1:ncol] - y$ES[(1+ncol):(2*ncol)])/2 * abs(sign(y$ES[1:ncol]) - sign(y$ES[(1+ncol):(2*ncol)]))/2)
      es <- t(matrix(unlist(zes), ncol = dim(ds2@mat)[2]))
      
      rownames(es) <- ds2@cid
      colnames(es) <- ds1@cid
      return(es)
    } else if (metric == "pearson") {
      
      return(calcSimBlock(ds1, ds2, metric="pearson", parallel=0))
      
    } else if (metric == "spearman") {
      
      return(calcSimBlock(ds1, ds2, metric="spearman", parallel=0))
      
    }
  }
}
