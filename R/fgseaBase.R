#' fgseaBase
#' 
#' This function is a modification of the fgsea function that dispenses with the adaptive permutation testing
#' and doesn't return or calculate leading edge genes. This makes it much faster when the only objective is 
#' calculation of the raw GSEA enrichment score (e.g. for use as a distance function).
#' @param pathways List of gene sets to check
#' @param stats Named vector of gene-level stats. Names should be the same as in 'pathways'
#' @param nperm Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param scoreType This parameter defines the GSEA score type. Possible options are ("std", "pos", "neg")
#' @param nproc If not equal to zero sets BPPARAM to use nproc workers (default = 0).
#' @param gseaParam GSEA parameter value, all gene-level statis are raised to the power of 'gseaParam' before calculation of GSEA enrichment scores.
#' @param BPPARAM Parallelization parameter used in bplapply. Can be used to specify cluster to run. If not initialized explicitly or by setting 'nproc' default value 'bpparam()' is used.
#'
#' @return A table with GSEA results. Each row corresponds to a tested pathway. The columns are:
#' \itemize{
#'   \item pathway -- name of the pathway
#'   \item ES -- enrichment score, the GSEA test statistic
#'   \item NES -- normalized enrichment score, the ES divided by the mean value for pathways of the same size. An NES of 1 [or -1] corresponds to the centre of the null distribution (for positively [or negatively] scored pathways)
#'   \item size -- number of genes in the pathway after removing genes not present in names(stats)
#' }
#'
#' @export
#' @import data.table
#' @import BiocParallel
#' @import fastmatch
#' @import stats
#' 


fgseaBase <- function(pathways,
                      stats,
                      nperm,
                      minSize   = 1,
                      maxSize   = Inf,
                      scoreType = c("std", "pos", "neg"),
                      nproc     = 0,
                      gseaParam = 1,
                      BPPARAM   = NULL) {
  scoreType <- match.arg(scoreType)
  pp <- fgsea:::preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, scoreType)
  pathwaysFiltered <- pp$filtered
  pathwaysSizes <- pp$sizes
  stats <- pp$stats
  
  
  m <- length(pathwaysFiltered)
  
  if (m == 0) {
    return(data.table(pathway=character(),
                      pval=numeric(),
                      padj=numeric(),
                      ES=numeric(),
                      NES=numeric(),
                      nMoreExtreme=numeric(),
                      size=integer(),
                      leadingEdge=list()))
  }
  
  
  granularity <- max(1000, ceiling(nperm / 128)) # not having more than 128 threads
  permPerProc <- rep(granularity, floor(nperm / granularity))
  if (nperm - sum(permPerProc) > 0) {
    permPerProc <- c(permPerProc, nperm - sum(permPerProc))
  }
  
  pval=nLeZero=nGeZero=leZeroMean=geZeroMean=nLeEs=nGeEs=NULL
  
  seeds <- sample.int(10^9, length(permPerProc))
  
  BPPARAM <- fgsea:::setUpBPPARAM(nproc=nproc, BPPARAM=BPPARAM)
  
  ES <- do.call(rbind,
                lapply(pathwaysFiltered, calcGseaStat,
                       stats=stats,
                       returnLeadingEdge=FALSE,
                       scoreType=scoreType))
  
  gtZeroMean <- ifelse(is.na(mean(ES[ES > 0])), 1, mean(ES[ES > 0]))
  ltZeroMean <- ifelse(is.na(mean(ES[ES < 0])), 1, mean(ES[ES < 0]))
  NES <- numeric(length(ES))
  NES <- ES / ifelse(ES > 0, gtZeroMean, ltZeroMean)
  

  retdf <- data.frame(pathway=names(pathwaysFiltered), ES=ES, NES=NES, size=pathwaysSizes)
  
  retdf
}