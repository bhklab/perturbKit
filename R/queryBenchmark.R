#' queryBenchmark - Benchmark the performance of L1000 queries
#'
#' @param ds The dataset - a  GCT object (e.g. a gctx file from cmapR). The assumed structure has columns are samples, rows are genes.
#' @param k Numeric vector of values to evaluate, typically 10, 100, 1000
#' @param parallel Boolean, whether to evaluate the multicore performance from the library parallel.
#' @param numCores numeric, the number of processing cores to assign to the task. Only used if parallel=1.
#'
#' @return
#' @export
#' @importFrom rbenchmark benchmark
#' @importFrom cmapR subset_gct
#' @importFrom parallel detectCores
queryBenchmark <- function(ds, k=c(10, 30, 100, 300), parallel=1, numCores=detectCores()){
  
  rpt <- c()
  metrics <- c("cosine", "wtcs", "pearson", "spearman")
  for (mymet in metrics){
    rpt[[mymet]] <- cbind(metric=mymet, mybenchmark(k, ds, myfunc=compute_sim_block, name="", args=c(metric=mymet, parallel=parallel), reps=10))
  }
 
  rpt <- do.call("rbind", rpt) 
  return(rpt)
}




mybenchmark <- function(n, ds, myfunc, name, args, reps=20){
  mybm <- c()
  for (ii in n){
    ds_slice <- cmapR::subset_gct(ds, cid=seq(ii))
    mybm <- rbind(mybm, rbenchmark::benchmark(do.call(myfunc, c(ds_slice, ds_slice, args)), replications=reps))
  }
  mybm$test <- paste(name, n)
  mybm
}
