#' Drug Matrix Kidney Dataset
#'
#' DEGs between drug (max dose) and control samples in DrugMatrix kidney data.
#' There are 57 drugs and DEGs are mouse genes.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of 21419 genes ranked by logFC*-log10(adj.p.val).
#'
#' @format Nested list of genesets.
#' @source \url{https://ntp.niehs.nih.gov/data/drugmatrix}
"drugmatrix.kidney"

#' Drug Matrix Liver Dataset
#'
#' DEGs between drug (max dose) and control samples in DrugMatrix liver data.
#' There are 57 drugs and  DEGs are mouse genes.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of 21419 genes ranked by logFC*-log10(adj.p.val).
#'
#' @format Nested list of genesets
#' @source \url{https://ntp.niehs.nih.gov/data/drugmatrix}
"drugmatrix.liver"

#' Sciplex A549 Dataset
#'
#' DEGs between drug (max dose) and control samples in SciPlex A549 cell line data.
#' There are 60 drugs and DEGs are EnsemblIDs.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of genes ranked by logFC*-log10(adj.p.val).
#' The length of `up_full` ranges from 2851 - 5782 based on the perturbation.
#'
#' @format Nested list of genesets
#' @source \url{https://figshare.com/articles/dataset/sciPlex_dataset/24681285?file=43381398}
"sciplex.a549"

#' Sciplex K562 Dataset
#'
#' DEGs between drug (max dose) and control samples in SciPlex K562 cell line data.
#' There are 60 drugs and DEGs are EnsemblIDs.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of genes ranked by logFC*-log10(adj.p.val).
#' The length of `up_full` ranges from 2778 - 3776 based on the perturbation.
#'
#' @format Nested list of genesets
#' @source \url{https://figshare.com/articles/dataset/sciPlex_dataset/24681285?file=43381398}
"sciplex.k562"

#' Sciplex MCF7 Dataset
#'
#' DEGs between drug (max dose) and control samples in SciPlex MCF7 cell line data.
#' There are 60 drugs and DEGs are EnsemblIDs.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of genes ranked by logFC*-log10(adj.p.val).
#' The length of `up_full` ranges from 5704 - 7259 based on the perturbation.
#'
#' @format Nested list of genesets
#' @source \url{https://figshare.com/articles/dataset/sciPlex_dataset/24681285?file=43381398}
"sciplex.mcf7"

#' Perturb-seq K562 Dataset
#'
#' DEGs between drug (max dose) and control samples in K562 cell line data.
#' There are 642 CRISPRi knockdowns and DEGs are HGNC symbols.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of 8561 genes ranked by logFC*-log10(adj.p.val).
#' The single cell dataset these genesets were derived from had 8561 genes.
#'
#' @format Nested list of genesets
#' @source \url{https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_processed_Perturb-seq_datasets/20029387}
"perturbseq.k562"

#' Perturb-seq RPE1 Dataset
#'
#' DEGs between drug (max dose) and control samples in RPE1 cell line data.
#' There are 642 CRISPRi knockdowns and DEGs are HGNC symbols.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of 8748 genes ranked by logFC*-log10(adj.p.val).
#' The single cell dataset these genesets were derived from had 8748 genes.
#'
#' @format Nested list of genesets
#' @source \url{https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_processed_Perturb-seq_datasets/20029387}
"perturbseq.rpe1"

#' Neurips 2023 B Cell signatures
#'
#' DEGs between drug and control samples (DMSO) in human pbmcs.
#' There are 135 shared drugs across the 4 cell types and DEGs are HGNC symbols.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of 6349 genes ranked by logFC*-log10(adj.p.val).
#' The pseudobulk dataset these genesets were derived from had 6349 genes.
#'
#' @format Nested list of genesets
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945}
"neurips2023.b"

#' Neurips 2023 NK Cell signatures
#'
#' DEGs between drug and control samples (DMSO) in human NK cells.
#' There are 135 shared drugs across the 4 cell types and DEGs are HGNC symbols.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of 6349 genes ranked by logFC*-log10(adj.p.val).
#' The pseudobulk dataset these genesets were derived from had 6349 genes.
#'
#' @format Nested list of genesets
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945}
"neurips2023.nk"

#' Neurips 2023 T Cell signatures
#'
#' DEGs between drug and control samples (DMSO) in human T cells.
#' There are 135 shared drugs across the 4 cell types and DEGs are HGNC symbols.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of 6349 genes ranked by logFC*-log10(adj.p.val).
#' The pseudobulk dataset these genesets were derived from had 6349 genes.
#'
#' @format Nested list of genesets
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945}
"neurips2023.t"

#' Neurips 2023 Myeloid Cell signatures
#'
#' DEGs between drug and control samples (DMSO) in human myeloid cells.
#' There are 135 shared drugs across the 4 cell types and DEGs are HGNC symbols.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of 6349 genes ranked by logFC*-log10(adj.p.val).
#' The pseudobulk dataset these genesets were derived from had 6349 genes.
#'
#' @format Nested list of genesets
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945}
"neurips2023.myeloid"
