## Transform seurat to loom file

## Dependencies
library(Seurat)
library(loomR)

print('Calling R')

## Loading seurat file
seurat_file <- snakemake@input[[1]]
seurat <- readRDS(seurat_file)
seurat <- FindVariableFeatures(object = seurat)
print(seurat_file)

## Set output path
loom <- snakemake@output[[1]]
print(loom)

## Cleaning annotations
for (i in 1:ncol(seurat@meta.data)) {
        if ( class(seurat@meta.data[,i]) %in% c('character', 'factor')){
                seurat@meta.data[ ,i] <- as.character(seurat@meta.data[,i])
                seurat@meta.data[ ,i][is.na(seurat@meta.data[ ,i])] <- 'NA'
        }
}

seurat.loom <- as.loom(x = seurat, 
                      filename = loom,
                      overwrite = TRUE)
seurat.loom$close_all()
