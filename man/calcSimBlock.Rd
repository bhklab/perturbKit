% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcSimBlock.R
\name{calcSimBlock}
\alias{calcSimBlock}
\title{calcSimBlock}
\usage{
calcSimBlock(
  ds1,
  ds2,
  metric = "cosine",
  kgenes = 0,
  gseaParam = 1,
  nperms = 1,
  parallel = 0,
  numCores = parallel::detectCores()
)
}
\arguments{
\item{ds1}{A gct instance, dataset 1 or a matrix}

\item{ds2}{A gct instance, dataset 2 (could be the same as dataset 1) or a matrix (must be same datatype as ds1)}

\item{metric}{A string, one of wtcs, cosine, spearman, pearson, fastwtcs}

\item{kgenes}{For those metrics that restrict the number of genes (wtcs, cosine), a parameter that determines how many genes are used to compute similarity. Default is 50 for GSEA, 0 (all) for cosine.}

\item{gseaParam}{The GSEA parameter that affects the weighting of the gene values in GSEA.  Default is 1.}

\item{nperms}{Numeric, indicates how many permutations are run by fgsea for significance calculations.  As p-values from GSEA are not used, this defaults to 1.}

\item{parallel}{Logical, whether to parallelize. Default = 0.}

\item{numCores}{numeric. If parallel, how many cores to use. Default = parallel::detectCores()}
}
\value{
matrix
}
\description{
This function takes as input two gct datasets and computes a similarity matrix between
columns of the first dataset and columns of the second.  It requires the gene space to be
common to both.  Note that some similarity metrics, like GSEA-based metrics, are not symmetric.
ds1 and ds2 may also be matrices of dimension N_features x N_signatures.
}
