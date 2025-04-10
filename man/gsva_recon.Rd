% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/baseline_methods.R
\name{gsva_recon}
\alias{gsva_recon}
\title{Reconstruct Gene Signatures Using GSVA Scores}
\usage{
gsva_recon(sce, sigs)
}
\arguments{
\item{sce}{A SingleCellExperiment object containing gene expression data.}

\item{sigs}{A list of gene signatures, where each element is a character vector of gene names.}
}
\value{
A list of reconstructed gene signatures, with the same structure as the input \code{sigs}.
}
\description{
This function reconstructs gene signatures based on their correlation with GSVA (Gene Set Variation Analysis) scores.
}
\details{
The function performs the following steps:
\enumerate{
\item Calculates GSVA scores for the input signatures.
\item Computes the correlation between gene expression and GSVA scores.
\item Ranks genes based on their correlation with each signature's GSVA score.
\item Selects the top-ranking genes to form new signatures of the same length as the original ones.
}
}
\note{
This function requires the GSVA, stats, and dplyr packages.
}
