% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{perturbseq.rpe1}
\alias{perturbseq.rpe1}
\title{Perturb-seq RPE1 Dataset}
\format{
Nested list of genesets
}
\source{
\url{https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_processed_Perturb-seq_datasets/20029387}
}
\usage{
perturbseq.rpe1
}
\description{
DEGs between drug (max dose) and control samples in RPE1 cell line data.
There are 642 CRISPRi knockdowns and DEGs are HGNC symbols.
For each geneset, the \code{up} set represents the top 100 DEGs
and the \code{up_full} represents the full ranked list of 8748 genes ranked by logFC*-log10(adj.p.val).
The single cell dataset these genesets were derived from had 8748 genes.
}
\keyword{datasets}
