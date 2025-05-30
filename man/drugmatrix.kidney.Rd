% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{drugmatrix.kidney}
\alias{drugmatrix.kidney}
\title{Drug Matrix Kidney Dataset}
\format{
Nested list of genesets.
}
\source{
\url{https://ntp.niehs.nih.gov/data/drugmatrix}
}
\usage{
drugmatrix.kidney
}
\description{
DEGs between drug (max dose) and control samples in DrugMatrix kidney data.
There are 57 drugs and DEGs are mouse genes.
For each geneset, the \code{up} set represents the top 100 DEGs
and the \code{up_full} represents the full ranked list of 21419 genes ranked by logFC*-log10(adj.p.val).
}
\keyword{datasets}
