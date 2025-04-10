% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{rank.var.eset}
\alias{rank.var.eset}
\title{Rank Genes in an ExpressionSet by Variability}
\usage{
rank.var.eset(eset, fn = mad, filter_zero = FALSE)
}
\arguments{
\item{eset}{An ExpressionSet object containing gene expression data.}

\item{fn}{A function to measure variability. Default is \code{mad} (median absolute deviation).}

\item{filter_zero}{A logical indicating whether to filter out genes with zero variance. Default is FALSE.}
}
\value{
A pvector object containing ranked gene names.
}
\description{
This function ranks genes in an ExpressionSet object based on their variability,
as measured by a specified function (default is median absolute deviation).
}
\details{
The function calculates the variability of each gene using the specified function (default: mad).
Genes are then sorted in descending order of variability. If \code{filter_zero} is TRUE,
genes with zero variance are removed before ranking.
}
