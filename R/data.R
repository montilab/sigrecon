#' Drug Matrix Kidney Dataset
#'
#' DEGs between drug (max dose) and control samples in DrugMatrix kidney data.
#' There are 39 drugs and DEGs are mouse genes.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of 15248 genes ranked by logFC*-log10(adj.p.val).
#'
#' @format Nested list of genesets.
#' @source \url{https://ntp.niehs.nih.gov/data/drugmatrix}
"drugmatrix.kidney"

#' Drug Matrix Liver Dataset
#'
#' DEGs between drug (max dose) and control samples in DrugMatrix liver data.
#' There are 39 drugs and DEGs are mouse genes.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of 15248 genes ranked by logFC*-log10(adj.p.val).
#'
#' @format Nested list of genesets
#' @source \url{https://ntp.niehs.nih.gov/data/drugmatrix}
"drugmatrix.liver"

#' Sciplex A549 Dataset
#'
#' DEGs between drug (max dose) and control samples in SciPlex A549 cell line data.
#' There are 23 drugs and DEGs are EnsemblIDs.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of genes ranked by logFC*-log10(adj.p.val).
#' The length of `up_full` ranges from 2789 - 15469 based on the perturbation.
#'
#' @format Nested list of genesets
#' @source \url{https://figshare.com/articles/dataset/sciPlex_dataset/24681285?file=43381398}
"sciplex.a549"

#' Sciplex K562 Dataset
#'
#' DEGs between drug (max dose) and control samples in SciPlex K562 cell line data.
#' There are 23 drugs and DEGs are EnsemblIDs.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of genes ranked by logFC*-log10(adj.p.val).
#' The length of `up_full` ranges from 2789 - 14703 based on the perturbation.
#'
#' @format Nested list of genesets
#' @source \url{https://figshare.com/articles/dataset/sciPlex_dataset/24681285?file=43381398}
"sciplex.k562"

#' Sciplex MCF7 Dataset
#'
#' DEGs between drug (max dose) and control samples in SciPlex MCF7 cell line data.
#' There are 23 drugs and DEGs are EnsemblIDs.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of genes ranked by logFC*-log10(adj.p.val).
#' The length of `up_full` ranges from 6797 - 20234 based on the perturbation.
#'
#' @format Nested list of genesets
#' @source \url{https://figshare.com/articles/dataset/sciPlex_dataset/24681285?file=43381398}
"sciplex.mcf7"

#' Perturb-seq K562 Dataset
#'
#' DEGs between drug (max dose) and control samples in K562 cell line data.
#' There are 1297 CRISPRi knockdowns and DEGs are HGNC symbols.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of genes ranked by logFC*-log10(adj.p.val).
#' The length of `up_full` ranges from 429 - 8563 based on the perturbation.
#'
#' @format Nested list of genesets
#' @source \url{https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_processed_Perturb-seq_datasets/20029387}
"perturbseq.k562"

#' Perturb-seq RPE1 Dataset
#'
#' DEGs between drug (max dose) and control samples in RPE1 cell line data.
#' There are 1297 CRISPRi knockdowns and DEGs are HGNC symbols.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of genes ranked by logFC*-log10(adj.p.val).
#' The length of `up_full` ranges from 438 - 8749 based on the perturbation.
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

#' GTEX Whole Blood Aging signatures (Discrete Bins)
#'
#' DEGs between 70-79 and 20-29 whole blood bulk samples from GTEX v10.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of 23131 genes ranked by logFC*-log10(adj.p.val).
#'
#' @format Nested list of genesets
#' @source \url{https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression}
"gtex.blood.aging.20.70"

#' GTEX Whole Brain Aging signatures (Discrete Bins)
#'
#' DEGs between 70-79 and 20-29 brain bulk samples from GTEX v10.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of 23131 genes ranked by logFC*-log10(adj.p.val).
#'
#' @format Nested list of genesets
#' @source \url{https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression}
"gtex.brain.aging.20.70"

#' GTEX Whole Blood Aging signatures (Continuous)
#'
#' DEGs between old and young whole blood bulk samples (age as a numerical covariate) from GTEX v10.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of 23131 genes ranked by logFC*-log10(adj.p.val).
#'
#' @format Nested list of genesets
#' @source \url{https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression}
"gtex.blood.aging"

#' GTEX Whole Brain Aging signatures (Continuous)
#'
#' DEGs between old and young whole brain bulk samples (age as a numerical covariate) from GTEX v10.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of 23131 genes ranked by logFC*-log10(adj.p.val).
#'
#' @format Nested list of genesets
#' @source \url{https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression}
"gtex.brain.aging"

#' GTEX Aging signatures (Continuous)
#'
#' DEGs between old and young tissue bulk samples (age as a numerical covariate) from GTEX v10.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of 23131 genes ranked by logFC*-log10(adj.p.val).
#'
#' @format Nested list of genesets
#' @source \url{https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression}
"gtex.aging"

#' Tahoe A498 Dataset
#'
#' DEGs between drug and control samples in Tahoe A498 cell line data.
#' There are 108 drugs and gene identifiers are a mixture of HGNC symbols and Ensembl IDs.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of genes ranked by logFC*-log10(adj.p.val).
#' The length of `up_full` ranges from 3982 - 41479 based on the perturbation.
#'
#' @format Nested list of genesets
#' @source \url{https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M}
"tahoe.a498"

#' Tahoe HCT15 Dataset
#'
#' DEGs between drug and control samples in Tahoe HCT15 cell line data.
#' There are 105 drugs and gene identifiers are a mixture of HGNC symbols and Ensembl IDs.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of genes ranked by logFC*-log10(adj.p.val).
#' The length of `up_full` ranges from 2393 - 36874 based on the perturbation.
#'
#' @format Nested list of genesets
#' @source \url{https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M}
"tahoe.hct15"

#' Tahoe HEC-1-A Dataset
#'
#' DEGs between drug and control samples in Tahoe HEC-1-A cell line data.
#' There are 109 drugs and gene identifiers are a mixture of HGNC symbols and Ensembl IDs.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of genes ranked by logFC*-log10(adj.p.val).
#' The length of `up_full` ranges from 3979 - 28808 based on the perturbation.
#'
#' @format Nested list of genesets
#' @source \url{https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M}
"tahoe.hec_1_a"

#' Tahoe LOVO Dataset
#'
#' DEGs between drug and control samples in Tahoe LOVO cell line data.
#' There are 109 drugs and gene identifiers are a mixture of HGNC symbols and Ensembl IDs.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of genes ranked by logFC*-log10(adj.p.val).
#' The length of `up_full` ranges from 3968 - 27303 based on the perturbation.
#'
#' @format Nested list of genesets
#' @source \url{https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M}
"tahoe.lovo"

#' Tahoe MIAPACA-2 Dataset
#'
#' DEGs between drug and control samples in Tahoe MIAPACA-2 cell line data.
#' There are 109 drugs and gene identifiers are a mixture of HGNC symbols and Ensembl IDs.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of genes ranked by logFC*-log10(adj.p.val).
#' The length of `up_full` ranges from 7620 - 39047 based on the perturbation.
#'
#' @format Nested list of genesets
#' @source \url{https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M}
"tahoe.miapaca_2"

#' Tahoe NCI-H23 Dataset
#'
#' DEGs between drug and control samples in Tahoe NCI-H23 cell line data.
#' There are 108 drugs and gene identifiers are a mixture of HGNC symbols and Ensembl IDs.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of genes ranked by logFC*-log10(adj.p.val).
#' The length of `up_full` ranges from 5448 - 39771 based on the perturbation.
#'
#' @format Nested list of genesets
#' @source \url{https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M}
"tahoe.nci_h23"

#' Tahoe PANC03.27 Dataset
#'
#' DEGs between drug and control samples in Tahoe PANC03.27 cell line data.
#' There are 106 drugs and gene identifiers are a mixture of HGNC symbols and Ensembl IDs.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of genes ranked by logFC*-log10(adj.p.val).
#' The length of `up_full` ranges from 3122 - 38165 based on the perturbation.
#'
#' @format Nested list of genesets
#' @source \url{https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M}
"tahoe.panc03.27"

#' Tahoe SNU-1 Dataset
#'
#' DEGs between drug and control samples in Tahoe SNU-1 cell line data.
#' There are 105 drugs and gene identifiers are a mixture of HGNC symbols and Ensembl IDs.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of genes ranked by logFC*-log10(adj.p.val).
#' The length of `up_full` ranges from 2393 - 36250 based on the perturbation.
#'
#' @format Nested list of genesets
#' @source \url{https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M}
"tahoe.snu_1"

#' Tahoe SNU-423 Dataset
#'
#' DEGs between drug and control samples in Tahoe SNU-423 cell line data.
#' There are 108 drugs and gene identifiers are a mixture of HGNC symbols and Ensembl IDs.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of genes ranked by logFC*-log10(adj.p.val).
#' The length of `up_full` ranges from 3148 - 39444 based on the perturbation.
#'
#' @format Nested list of genesets
#' @source \url{https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M}
"tahoe.snu_423"

#' Tahoe SW48 Dataset
#'
#' DEGs between drug and control samples in Tahoe SW48 cell line data.
#' There are 104 drugs and gene identifiers are a mixture of HGNC symbols and Ensembl IDs.
#' For each geneset, the `up` set represents the top 100 DEGs
#' and the `up_full` represents the full ranked list of genes ranked by logFC*-log10(adj.p.val).
#' The length of `up_full` ranges from 2393 - 36402 based on the perturbation.
#'
#' @format Nested list of genesets
#' @source \url{https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M}
"tahoe.sw48"
