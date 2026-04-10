library(ArchR)

# cache arrow files from hise
library(hise)
files <- cacheFiles(list("04323ec7-7768-4617-80e1-9ad6cc191951", "f07e7b20-7e9c-440b-aaa3-82433faab11f", "e4b23f20-7283-4c35-b5bd-d01561581fa9", "7517621d-7505-4008-b818-0c4375723b60", "f8915652-9b0f-4501-9752-0293f6d78a4d", "375e21ab-154e-44dd-8fdf-e4a50ce4bbc0", "54ece30f-4576-47d3-8678-e227ba81ce09", "da3fae63-0fec-4083-8c8c-6b3e2bc5d720", "02b36580-3676-4559-96b9-1fae81167bde", "c1d47ab1-fc71-453c-9d14-0c929e191b67", "089ed0ea-96e6-4851-b204-757cfc574800", "eac97bc4-d52d-4677-8e7d-40fd2554fcb9"))

## Pull in just 6 files. The other 6 are duplicates. 
allArrows = list.files(path='input', pattern ='.arrow', recursive=TRUE, full.names=TRUE)
file_names <- basename(allArrows)
unique_indices <- !duplicated(file_names)
file_paths_unique <- allArrows[unique_indices]
file_names_unique <- file_names[unique_indices]
# names(file_paths_unique) <- file_names_unique

addArchRGenome('hg38')
proj <- ArchRProject(file_paths_unique, outputDirectory = 'Plasma_ArchR')
saveArchRProject(proj)

proj <- loadArchRProject('Plasma_ArchR')
metadf <- as.data.frame(getCellColData(proj))

### Load Seurat object metadata. 
seur <- readRDS('teaseq_4_05.rds')
metadt <- as.data.frame(seur@meta.data)

table(metadt$meta)
table(metadt$subisotype)
table(metadt$isotype)
table(metadt$isotype, metadt$pbmc_sample_id)
table(metadt$subisotype, metadt$pbmc_sample_id)
## There's enough cells per sample to call peaks on sub-isotypes
sum(!metadt$barcodes %in% metadf$barcodes)
# 2597 cells in seurat object not in arrow files

head(metadt[,c('barcodes', 'pbmc_sample_id', 'pool_id', 'batch_id')])
metadt$atac_barcodes = paste(metadt$pool_id, "_", metadt$pbmc_sample_id, "#", metadt$barcodes, sep = '')
sum(metadt$atac_barcodes %in% rownames(metadf))/dim(metadt)[1]
table(!metadt$atac_barcodes %in% rownames(metadf), metadt$pbmc_sample_id)

metadf$Isotype = metadt$isotype[match(rownames(metadf), 
                                        metadt$atac_barcodes)]
metadf$Subisotype = metadt$subisotype[match(rownames(metadf), 
                                        metadt$atac_barcodes)]
metadf$Meta = metadt$meta[match(rownames(metadf), 
                                        metadt$atac_barcodes)]
proj <- addCellColData(proj, data = metadf$Isotype, name = 'Isotype', 
 cells = rownames(metadf))
proj <- addCellColData(proj, data = metadf$Subisotype, name = 'Subisotype', 
    cells = rownames(metadf))
proj <- addCellColData(proj, data = metadf$Meta, name = 'Meta', 
    cells = rownames(metadf))
saveArchRProject(proj)
