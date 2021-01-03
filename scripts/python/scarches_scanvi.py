## Scarches SCANVI implementation

## Dependencies
import scarches as sca
import scanpy as sc
import pandas as pd
import anndata as ad
import os
import warnings
import torch
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt
import numpy as np
import gdown
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
sc.settings.set_figure_params(dpi=200, frameon=False)
sc.set_figure_params(dpi=200)
sc.set_figure_params(figsize=(4, 4))
torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7) 


##############################################################
## Definition of Auxilliary functions

## Subsetting adata by top HVGs
def subset_hvg(adata, n_top):
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top)
    adata = adata[:, adata.var['highly_variable']]
    return(adata)


#############################################################

## Reading reference and query datasets
annData = sc.read_loom(snakemake.input.reference_loom)
target_annData = sc.read_loom(snakemake.input.query_loom)


#############################################################
##							   ##
## 	Harmonizing reference and query datasets	   ##
##							   ##
#############################################################  

## Subsetting annotations to cell type and batches
annData.obs = annData.obs[
                 [snakemake.params.ref_cell_type_column, 
                 snakemake.params.ref_batch_column]
]

## Harmonizing query with reference annotations
renamed_df = target_annData.obs.rename(
     columns={snakemake.params.query_batch_column:snakemake.params.ref_batch_column, 
              snakemake.params.query_cell_type_column:snakemake.params.ref_cell_type_column}
)
target_annData.obs = renamed_df

## Subsetting to top 5000 HGVs
annData = subset_hvg(annData, snakemake.params.n_top_hvg)
target_annData = subset_hvg(target_annData, snakemake.params.n_top_hvg)

## Subsampling
sc.pp.subsample(annData, snakemake.params.subsampling_fraction)

## Subset to genes in the intersection
genes_int = set(annData.var_names).intersection(set(target_annData.var_names))
genes_int = list(genes_int)
reference = annData[:, genes_int].copy()
query = target_annData[:, genes_int].copy()


##############################################################################
##									    ##
##			Running SCANVI					    ##
##									    ##
##############################################################################

## NN settings Settings
vae_epochs = 20
scanvi_epochs = 20
surgery_epochs = 20

early_stopping_kwargs = {
    "early_stopping_metric": "elbo",
    "save_best_state_metric": "elbo",
    "patience": 10,
    "threshold": 0,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}
early_stopping_kwargs_scanvi = {
    "early_stopping_metric": "accuracy",
    "save_best_state_metric": "accuracy",
    "on": "full_dataset",
    "patience": 10,
    "threshold": 0.001,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}
early_stopping_kwargs_surgery = {
    "early_stopping_metric": "elbo",
    "save_best_state_metric": "elbo",
    "on": "full_dataset",
    "patience": 10,
    "threshold": 0.001,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}

## Setting up data
sca.dataset.setup_anndata(reference, 
                          batch_key=snakemake.params.ref_batch_column, 
                          labels_key=snakemake.params.ref_cell_type_column)

## Definition of the model
vae = sca.models.SCANVI(
    reference,
    "Unknown",
    n_layers=2,
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
)

## Trainning the model
vae.train(
    n_epochs_unsupervised=snakemake.params.vae_epochs,
    n_epochs_semisupervised=snakemake.params.scanvi_epochs,
    unsupervised_trainer_kwargs=dict(early_stopping_kwargs=early_stopping_kwargs),
    semisupervised_trainer_kwargs=dict(metrics_to_monitor=["elbo", "accuracy"],
                                       early_stopping_kwargs=early_stopping_kwargs_scanvi),
    frequency=1,
    n_epochs_kl_warmup=snakemake.params.n_epochs_kl_warmup,
)
print('Training finished')

## Calculate latent representation
reference_latent = sc.AnnData(vae.get_latent_representation())
reference_latent.obs["cell_type"] = reference.obs[snakemake.params.ref_cell_type_column].tolist()
reference_latent.obs["batch"] = reference.obs[snakemake.params.ref_batch_column].tolist()

## Projecting latent representation into UMAP
sc.pp.neighbors(reference_latent, n_neighbors=8)
sc.tl.leiden(reference_latent)
sc.tl.umap(reference_latent)

## UMAP plot
umap = sc.pl.umap(reference_latent,
           color=['batch', 'cell_type'],
           frameon=False,
           wspace=0.6,
           save=snakemake.params.reference_latent_umap,
           )

## Calculating accuracy
reference_latent.obs['predictions'] = vae.predict()



##################################################################################
##										##
##		Surgery on reference model and train on query			##
##										##
##################################################################################

print(query.obs['Cluster'][0:2])

model = sca.models.SCANVI.load_query_data(
    query,
    vae,
    freeze_dropout = True,
)
model._unlabeled_indices = np.arange(query.n_obs)
model._labeled_indices = []
print("Labelled Indices: ", len(model._labeled_indices))
print("Unlabelled Indices: ", len(model._unlabeled_indices))


## Training the integrated model
model.train(
    n_epochs_semisupervised=surgery_epochs,
    train_base_model=False,
    semisupervised_trainer_kwargs=dict(metrics_to_monitor=["accuracy", "elbo"],
                                       weight_decay=0,
                                       early_stopping_kwargs=early_stopping_kwargs_surgery
                                      ),
    frequency=1
)

## Getting latent space
query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['predictions'] = model.predict()
sc.pp.neighbors(query_latent)
sc.tl.leiden(query_latent)
sc.tl.umap(query_latent)
sc.pl.umap(
    query_latent,
    color=["predictions"],
    frameon=False,
    wspace=0.6,
    save=snakemake.params.query_latent_umap,
)

## Getting latent representation of integrated data
adata_full = reference.concatenate(query)
full_latent = sc.AnnData(model.get_latent_representation(adata=adata_full))
full_latent.obs['cell_type'] = adata_full.obs[snakemake.params.ref_cell_type_column].tolist()
full_latent.obs['batch'] = adata_full.obs[snakemake.params.ref_batch_column].tolist()
sc.pp.neighbors(full_latent)
sc.tl.leiden(full_latent)
sc.tl.umap(full_latent)
plt.figure()
sc.pl.umap(
    full_latent,
    color=["batch", "cell_type"],
    frameon=False,
    wspace=0.6,
    save=snakemake.params.integrated_latent_umap,
)

print('Done')
