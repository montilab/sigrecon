# List available on-demand datasets

Returns a table of the perturbational signature datasets available via
[`get_dataset()`](https://montilab.github.io/sigrecon/reference/get_dataset.md),
with a short description and original data source for each.

## Usage

``` r
list_datasets()
```

## Value

A data.frame with columns `name`, `description`, `source`.

## Examples

``` r
list_datasets()
#>                   name
#> 1    drugmatrix.kidney
#> 2     drugmatrix.liver
#> 3         sciplex.a549
#> 4         sciplex.k562
#> 5         sciplex.mcf7
#> 6      perturbseq.k562
#> 7      perturbseq.rpe1
#> 8        neurips2023.b
#> 9       neurips2023.nk
#> 10       neurips2023.t
#> 11 neurips2023.myeloid
#> 12          gtex.blood
#> 13    gtex.brain.hippo
#> 14          tahoe.a498
#> 15         tahoe.hct15
#> 16       tahoe.hec_1_a
#> 17          tahoe.lovo
#> 18     tahoe.miapaca_2
#> 19       tahoe.nci_h23
#> 20     tahoe.panc03.27
#> 21         tahoe.snu_1
#> 22       tahoe.snu_423
#> 23          tahoe.sw48
#>                                                                                                                    description
#> 1                           DEGs between drug (max dose) and control samples in DrugMatrix kidney data. 39 drugs, mouse genes.
#> 2                            DEGs between drug (max dose) and control samples in DrugMatrix liver data. 39 drugs, mouse genes.
#> 3                       DEGs between drug (max dose) and control samples in SciPlex A549 cell line data. 23 drugs, EnsemblIDs.
#> 4                       DEGs between drug (max dose) and control samples in SciPlex K562 cell line data. 23 drugs, EnsemblIDs.
#> 5                       DEGs between drug (max dose) and control samples in SciPlex MCF7 cell line data. 23 drugs, EnsemblIDs.
#> 6  DEGs between drug (max dose) and control samples in Perturb-seq K562 cell line data. 1297 CRISPRi knockdowns, HGNC symbols.
#> 7  DEGs between drug (max dose) and control samples in Perturb-seq RPE1 cell line data. 1297 CRISPRi knockdowns, HGNC symbols.
#> 8             DEGs between drug and control (DMSO) samples in NeurIPS 2023 human PBMC B cells. 135 shared drugs, HGNC symbols.
#> 9            DEGs between drug and control (DMSO) samples in NeurIPS 2023 human PBMC NK cells. 135 shared drugs, HGNC symbols.
#> 10            DEGs between drug and control (DMSO) samples in NeurIPS 2023 human PBMC T cells. 135 shared drugs, HGNC symbols.
#> 11      DEGs between drug and control (DMSO) samples in NeurIPS 2023 human PBMC myeloid cells. 135 shared drugs, HGNC symbols.
#> 12                             DEGs between old and young whole blood bulk samples (age as numerical covariate) from GTEx v10.
#> 13                       DEGs between old and young brain hippocampus bulk samples (age as numerical covariate) from GTEx v10.
#> 14                            DEGs between drug and control samples in Tahoe A498 cell line data. 108 drugs, HGNC/Ensembl mix.
#> 15                           DEGs between drug and control samples in Tahoe HCT15 cell line data. 105 drugs, HGNC/Ensembl mix.
#> 16                         DEGs between drug and control samples in Tahoe HEC-1-A cell line data. 109 drugs, HGNC/Ensembl mix.
#> 17                            DEGs between drug and control samples in Tahoe LOVO cell line data. 109 drugs, HGNC/Ensembl mix.
#> 18                       DEGs between drug and control samples in Tahoe MIAPACA-2 cell line data. 109 drugs, HGNC/Ensembl mix.
#> 19                         DEGs between drug and control samples in Tahoe NCI-H23 cell line data. 108 drugs, HGNC/Ensembl mix.
#> 20                       DEGs between drug and control samples in Tahoe PANC03.27 cell line data. 106 drugs, HGNC/Ensembl mix.
#> 21                           DEGs between drug and control samples in Tahoe SNU-1 cell line data. 105 drugs, HGNC/Ensembl mix.
#> 22                         DEGs between drug and control samples in Tahoe SNU-423 cell line data. 108 drugs, HGNC/Ensembl mix.
#> 23                            DEGs between drug and control samples in Tahoe SW48 cell line data. 104 drugs, HGNC/Ensembl mix.
#>                                                                                                                                                                                          source
#> 1                                                                                                                                                     https://ntp.niehs.nih.gov/data/drugmatrix
#> 2                                                                                                                                                     https://ntp.niehs.nih.gov/data/drugmatrix
#> 3                                                                                                                  https://figshare.com/articles/dataset/sciPlex_dataset/24681285?file=43381398
#> 4                                                                                                                  https://figshare.com/articles/dataset/sciPlex_dataset/24681285?file=43381398
#> 5                                                                                                                  https://figshare.com/articles/dataset/sciPlex_dataset/24681285?file=43381398
#> 6  https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_processed_Perturb-seq_datasets/20029387
#> 7  https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_processed_Perturb-seq_datasets/20029387
#> 8                                                                                                                                  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945
#> 9                                                                                                                                  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945
#> 10                                                                                                                                 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945
#> 11                                                                                                                                 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945
#> 12                                                                                                                      https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
#> 13                                                                                                                      https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
#> 14                                                                                                                  https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M
#> 15                                                                                                                  https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M
#> 16                                                                                                                  https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M
#> 17                                                                                                                  https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M
#> 18                                                                                                                  https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M
#> 19                                                                                                                  https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M
#> 20                                                                                                                  https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M
#> 21                                                                                                                  https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M
#> 22                                                                                                                  https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M
#> 23                                                                                                                  https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M
```
