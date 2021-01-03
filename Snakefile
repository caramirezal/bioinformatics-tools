## Snakemake pipeline containing routines for high throughput sequencing 
## including RNA, ATAC-Seq bulk and single cell as well as multiomics.
configfile: 'config/config.yml' 


## Convert seurat formats to loom to load into python based algorithms
rule seurat2loom:
    input:
        config['seurat2loom']['seurat_file']
    output:
        config['seurat2loom']['out_loom_path']
    script:
        'scripts/R/seurat2loom.R'


## SCANVI scarches 
rule scarches_scanvi:
    input:
        reference_loom=config['scarches_scanvi']['reference_loom'],
        query_loom=config['scarches_scanvi']['query_loom']
    params:
        ref_batch_column=config['scarches_scanvi']['ref_batch_column'],
        ref_cell_type_column=config['scarches_scanvi']['ref_cell_type_column'],
        query_batch_column=config['scarches_scanvi']['query_batch_column'],
        query_cell_type_column=config['scarches_scanvi']['query_cell_type_column'],
        n_top_hvg=config['scarches_scanvi']['n_top_hvg'],
        subsampling_fraction=config['scarches_scanvi']['subsampling_fraction'],
        vae_epochs=config['scarches_scanvi']['vae_epochs'],
        scanvi_epochs=config['scarches_scanvi']['scanvi_epochs'],
        surgery_epochs=config['scarches_scanvi']['surgery_epochs'],
        n_epochs_kl_warmup=config['scarches_scanvi']['n_epochs_kl_warmup'],
        reference_latent_umap= '_' + config['scarches_scanvi']['reference_latent_umap'],
        integrated_latent_umap= '_' + config['scarches_scanvi']['integrated_latent_umap'],
        query_latent_umap= '_' + config['scarches_scanvi']['query_latent_umap']
    script:
        'scripts/python/scarches_scanvi.py'


