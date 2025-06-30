import scanpy as sc
import anndata as ad
from adpbulk import ADPBulk
from tqdm import tqdm
import gc

# Download raw Tahoe data from google cloud storage 
# See this for reference (https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M)
# On SCC /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/perturbational_data/tahoe
for i in tqdm(range(4, 14)):
    id = str(i + 1)
    PATH = f"plate{id}_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad"
    adata = ad.read_h5ad(PATH)
    adata.raw = adata

    adpb = ADPBulk(adata, ["sample", "cell_name", "drugname_drugconc"], use_raw=True)
    pseudobulk_matrix = adpb.fit_transform()
    pseudo_meta = adpb.get_meta()

    pseudo_adata = ad.AnnData(
        X=pseudobulk_matrix.values,
        obs=pseudo_meta,
        var=adata.var.loc[pseudobulk_matrix.columns].copy()
    )

    output_path = f"plate{id}_pseudobulk.h5ad"
    pseudo_adata.write_h5ad(output_path)

    del adata
    del pseudo_adata
    gc.collect()


data_dict = {}

for i in range(14):
    data_dict.update({f"plate_{i+1}": f"plate{i+1}_pseudobulk.h5ad"})

ad.experimental.concat_on_disk(
    data_dict,
    f"merged_pseudobulk.h5ad",
    label='plate',
)


# Other changes
# evaluate_dup_counts(ad_pb)
# There are no duplicated counts.
ad_pb = ad.read_h5ad("merged_pseudobulk.h5ad")
ad_pb.X = sparse.csr_matrix(ad_pb.X)
ad_pb.obs["drugname_drugconc"] = ad_pb.obs["drugname_drugconc"].apply(lambda x: x.replace("[","").replace("]","").replace("(","").replace(")","").replace("'",""))
ad_pb.obs["drug_name"] = ad_pb.obs["drugname_drugconc"].apply(lambda x: x.split(",")[0])
ad_pb.obs["drug_conc"] = ad_pb.obs["drugname_drugconc"].apply(lambda x: x.split(",")[1] + x.split(",")[2])
ad_pb.obs_names = ad_pb.obs['plate'].astype(str) + '_' + ad_pb.obs['cell_name'].astype(str) + "_" + ad_pb.obs['drug_name'].astype(str) + "_" + ad_pb.obs['drug_conc'].astype(str) + "_" + ad_pb.obs_names.astype(str)
for col in ad_pb.obs.select_dtypes(['category']).columns:
    ad_pb.obs[col] = ad_pb.obs[col].astype(str)
for col in ad_pb.var.select_dtypes(['category']).columns:
    ad_pb.var[col] = ad_pb.var[col].astype(str)
ad_pb.write_h5ad("merged_pseudobulk.h5ad")
