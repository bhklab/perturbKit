# Compute distance is a function for computing similarities 
compute_distance <- function(ds, gene_annot, upgenes=c(), dngenes=c(), geneweights=0, 
                             geneset_name="", metric="cmap_score", 
                             gseaParam=1, nperms=1000){
  
  # Add condition to check input, upgenes, dngenes, or a weighted vector geneweights
  # Add condition to check metric, currently only cmap_score supported. 
  
  res <- data.frame(cid=character(length(ds@cid)), 
                    geneset=character(length(ds@cid)), 
                    metric=character(length(ds@cid)), 
                    value=numeric(length(ds@cid)), 
                    value_norm=numeric(length(ds@cid)),
                    pval=numeric(length(ds@cid)), 
                    padj=numeric(length(ds@cid)), 
                    size=numeric(length(ds@cid)), stringsAsFactors=FALSE)
  
  if (metric == "cmap_score") {
    upgenes <- intersect(upgenes, gene_annot$pr_gene_symbol)
    dngenes <- intersect(dngenes, gene_annot$pr_gene_symbol)
    common <- intersect(upgenes, dngenes)
    
    upgenes <- setdiff(upgenes, common)
    dngenes <- setdiff(dngenes, common)
    
    upid <- geneinfo$pr_gene_id[match(upgenes, geneinfo$pr_gene_symbol)]
    dnid <- geneinfo$pr_gene_id[match(dngenes, geneinfo$pr_gene_symbol)]
    
    for (ii in seq(length(ds@cid))){
      if (ii %% 1000 == 0) {print(ii)}
      #tGSEA <- fgsea(list(upid=upid, dnid=dnid), ds@mat[, ii], nperm=nperms, gseaParam=gseaParam)
      tGSEA <- fgseaMultilevel(list(upid=upid, dnid=dnid), ds@mat[, ii], sampleSize=k, gseaParam) 
      combES <- (tGSEA$ES[1] - tGSEA$ES[2])/2 * abs(sign(tGSEA$ES[1]) - sign(tGSEA$ES[2]))/2
      combNES <- (tGSEA$NES[1] - tGSEA$NES[2])/2 * abs(sign(tGSEA$NES[1]) - sign(tGSEA$NES[2]))/2
      
      res[ii,] <- list(cid=ds@cid[ii], 
                       geneset=geneset_name, 
                       metric=metric,
                       value=combES, 
                       value_norm=combNES, 
                       pval=1,   #Fix, probably Fisher's method from fgseaMultilevel
                       padj=1,
                       size=tGSEA$size[1] + tGSEA$size[2])
    }
  }
  
  else if (metric == "cosine"){
    # 
    
  }
  
  else if (metric == "pearson"){
    
  }
  
  else if (metric == "spearman"){
    
  }
  
  else if (metric == "chardir"){
    
  }
  
  else if (metric == "wcosine"){
    
  }
  
  res
}