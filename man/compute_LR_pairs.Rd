% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_LR_pairs.R
\name{compute_LR_pairs}
\alias{compute_LR_pairs}
\title{Compute ligand-receptor pairs}
\usage{
compute_LR_pairs(
  RNA.tpm,
  remove.genes.ICB_proxies = FALSE,
  compute.cytokines.pairs = FALSE,
  cancertype
)
}
\arguments{
\item{RNA.tpm}{numeric matrix of tpm values with rows=genes and columns=samples}

\item{remove.genes.ICB_proxies}{boolean variable to reomove all those genes involved in the computation of ICB proxy's of response}

\item{compute.cytokines.pairs}{boolean variable to compute cytokine pairs as well}

\item{cancerype}{string character}
}
\value{
A list with the following elements:
\describe{
\item{LRpairs}{Ligand-leceptor pairs weights matrix in log2(tpm + 1) with rows=samples and columns = L-R pairs}
\item{CYTOKINEpairs}{Cytokine-Cytokine pairs weights matrix in log2(tpm + 1) with rows=samples and columns = CYTOKINE pairs}
}
}
\description{
\code{compute_LR_pairs} obtain ligand-receptor pairs weights from tpm RNAseq data.
}
