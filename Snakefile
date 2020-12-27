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



