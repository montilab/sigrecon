# Demo datasets

`sigrecon` bundles four small, real-data recontextualization demos (one
per major dataset covered by `scripts/`), each a `demo_<dataset>_se` /
`demo_<dataset>_sigs` / `demo_<dataset>_true_sigs` triple with the same
shape: `demo_<dataset>_se` is target-context expression,
`demo_<dataset>_sigs` are real source-context signatures (as if defined
in a different biological context), and `demo_<dataset>_true_sigs` is
the real target-context ground truth, for benchmarking
recontextualization quality with
[`sig_eval_table()`](https://montilab.github.io/sigrecon/reference/sig_eval_table.md).
Gene sets are restricted to the union of genes needed by the kept
signatures (not variance-based selection – see `data-raw/demo_sciplex.R`
for why), and perturbations are sampled down to ~23 for speed. All are
fast, self-contained examples for
[`projectCor()`](https://montilab.github.io/sigrecon/reference/projectCor.md),
[`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)/[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md),
and
[`sig_eval_table()`](https://montilab.github.io/sigrecon/reference/sig_eval_table.md)
– see the package README.

No data leakage: each `demo_<dataset>_se`'s perturbation samples are a
set of perturbations DISJOINT from the ~23 evaluated in the
corresponding `_sigs`/`_true_sigs` (control samples plus a *different*
set of ~23 "background" perturbations, not the evaluated ones) – so the
network/projection input never contains a sample for the perturbation
being predicted. See the "No data leakage" note in each
`data-raw/demo_*.R` script.

|  |  |  |  |
|----|----|----|----|
| Dataset | Source context | Target context (`_se`) | Built by |
| SciPlex | K562 | A549 | `data-raw/demo_sciplex.R` |
| Tahoe | NCI-H23 | A498 | `data-raw/demo_tahoe.R` |
| DrugMatrix | kidney | liver | `data-raw/demo_drugmatrix.R` |
| Perturb-seq | K562 | RPE1 | `data-raw/demo_perturbseq.R` |

## Usage

``` r
demo_sciplex_se

demo_sciplex_sigs

demo_sciplex_true_sigs

demo_tahoe_se

demo_tahoe_sigs

demo_tahoe_true_sigs

demo_drugmatrix_se

demo_drugmatrix_sigs

demo_drugmatrix_true_sigs

demo_perturbseq_se

demo_perturbseq_sigs

demo_perturbseq_true_sigs
```

## Format

`demo_sciplex_se`: a `SummarizedExperiment` with 1025 genes (Ensembl
IDs) and 48 samples (Vehicle control + 23 background drugs, disjoint
from the 23 evaluated in `demo_sciplex_sigs`/ `demo_sciplex_true_sigs`),
built from real SciPlex A549 pseudobulk expression
(<https://figshare.com/articles/dataset/sciPlex_dataset/24681285?file=43381398>).
Has `logcounts` (log2 CPM, default assay) and `counts` (raw) assays.
`colData` has `orig.ident`, `product_name` (drug, or "Vehicle"),
`replicate`.

`demo_sciplex_sigs`: named list (by drug) of character vectors of
Ensembl gene IDs – subset of `sciplex.k562`.

`demo_sciplex_true_sigs`: named list (by drug), each element a list with
`up`/`up_full` – subset of `sciplex.a549`.

`demo_tahoe_se`: a `SummarizedExperiment` with 1808 genes (HGNC symbols)
and 55 samples (DMSO_TF control + 23 background drugs, disjoint from the
23 evaluated in `demo_tahoe_sigs`/ `demo_tahoe_true_sigs`), built from
real Tahoe A498 pseudobulk expression
(<https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M>).
Has `logcounts` (log2 CPM, default assay) and `counts` (raw) assays.

`demo_tahoe_sigs`: named list (by drug) of character vectors – subset of
`tahoe.nci_h23` (source context).

`demo_tahoe_true_sigs`: named list (by drug), each element a list with
`up`/`up_full` – subset of `tahoe.a498` (target context).

`demo_drugmatrix_se`: a `SummarizedExperiment` with 1846 genes (gene
symbols) and 84 samples (15 control + max-dose/time samples for 23
background drugs, disjoint from the 23 evaluated in
`demo_drugmatrix_sigs`/`demo_drugmatrix_true_sigs`), built from real
DrugMatrix liver microarray expression
(<https://ntp.niehs.nih.gov/data/drugmatrix>). Single `logcounts` assay
(already log-transformed microarray intensities, no further
normalization applied). `colData` has `compound`, `dose`, `tissue`,
`time`, `vehicle`.

`demo_drugmatrix_sigs`: named list (by drug) of character vectors –
subset of `drugmatrix.kidney` (source context).

`demo_drugmatrix_true_sigs`: named list (by drug), each element a list
with `up`/`up_full` – subset of `drugmatrix.liver` (target context).

`demo_perturbseq_se`: a `SummarizedExperiment` with 1382 genes (HGNC
symbols) and 84 samples (15 non-targeting control + up to 3 replicates
each for 23 background CRISPRi knockdowns, disjoint from the 23
evaluated in `demo_perturbseq_sigs`/ `demo_perturbseq_true_sigs`), built
from real Perturb-seq RPE1 pseudobulk expression
(<https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_processed_Perturb-seq_datasets/20029387>).
Has `logcounts` (log2 CPM, default assay) and `counts` (raw) assays.
`colData` has `gene` (knockdown target, or "non-targeting"),
`gem_group`, `is_control`.

`demo_perturbseq_sigs`: named list (by knockdown target) of character
vectors – subset of `perturbseq.k562` (source context).

`demo_perturbseq_true_sigs`: named list (by knockdown target), each
element a list with `up`/`up_full` – subset of `perturbseq.rpe1` (target
context).
