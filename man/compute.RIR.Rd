% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute.RIR.R
\name{compute.RIR}
\alias{compute.RIR}
\title{Compute Immune resistance program}
\usage{
compute.RIR(RNA.tpm)
}
\arguments{
\item{RNA.tpm}{numeric matrix with rows=genes and columns=samples}
}
\value{
numeric matrix with rows=samples and columns=IRP
}
\description{
\code{compute.RIR} computes immune resistance program using the code provided by the authors.
(Jerby-Arnon et al., 2018).
}
