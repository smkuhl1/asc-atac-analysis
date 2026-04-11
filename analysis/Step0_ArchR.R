## Step 0: Initialize ArchR project and attach Seurat-derived metadata.
##
## Inputs:
## - Arrow files under input/
## - Seurat object: teaseq_4_05.rds
##
## Output:
## - ArchR project directory: Plasma_ArchR/

library(ArchR)
library(hise)

project_dir <- "Plasma_ArchR"
seurat_rds <- "teaseq_4_05.rds"

# Cache Arrow files from HiSE.
arrow_ids <- c(
    "04323ec7-7768-4617-80e1-9ad6cc191951", "f07e7b20-7e9c-440b-aaa3-82433faab11f",
    "e4b23f20-7283-4c35-b5bd-d01561581fa9", "7517621d-7505-4008-b818-0c4375723b60",
    "f8915652-9b0f-4501-9752-0293f6d78a4d", "375e21ab-154e-44dd-8fdf-e4a50ce4bbc0",
    "54ece30f-4576-47d3-8678-e227ba81ce09", "da3fae63-0fec-4083-8c8c-6b3e2bc5d720",
    "02b36580-3676-4559-96b9-1fae81167bde", "c1d47ab1-fc71-453c-9d14-0c929e191b67",
    "089ed0ea-96e6-4851-b204-757cfc574800", "eac97bc4-d52d-4677-8e7d-40fd2554fcb9"
)
files <- cacheFiles(as.list(arrow_ids))

# Keep one Arrow path per unique basename because half are duplicates.
all_arrows <- list.files(path = "input", pattern = "\\.arrow$", recursive = TRUE, full.names = TRUE)
arrow_names <- basename(all_arrows)
unique_arrows <- all_arrows[!duplicated(arrow_names)]

addArchRGenome("hg38")
proj <- ArchRProject(unique_arrows, outputDirectory = project_dir)
saveArchRProject(proj)

proj <- loadArchRProject(project_dir)
archr_meta <- as.data.frame(getCellColData(proj))

# Load Seurat metadata.
seur <- readRDS(seurat_rds)
seurat_meta <- as.data.frame(seur@meta.data)

table(seurat_meta$meta)
table(seurat_meta$subisotype)
table(seurat_meta$isotype)
table(seurat_meta$isotype, seurat_meta$pbmc_sample_id)
table(seurat_meta$subisotype, seurat_meta$pbmc_sample_id)

# Expected mismatch check between Seurat and Arrow barcodes.
sum(!seurat_meta$barcodes %in% archr_meta$barcodes)

head(seurat_meta[, c("barcodes", "pbmc_sample_id", "pool_id", "batch_id")])
seurat_meta$atac_barcodes <- paste(
    seurat_meta$pool_id,
    "_",
    seurat_meta$pbmc_sample_id,
    "#",
    seurat_meta$barcodes,
    sep = ""
)

sum(seurat_meta$atac_barcodes %in% rownames(archr_meta)) / nrow(seurat_meta)
table(!seurat_meta$atac_barcodes %in% rownames(archr_meta), seurat_meta$pbmc_sample_id)

archr_meta$Isotype <- seurat_meta$isotype[match(rownames(archr_meta), seurat_meta$atac_barcodes)]
archr_meta$Subisotype <- seurat_meta$subisotype[match(rownames(archr_meta), seurat_meta$atac_barcodes)]
archr_meta$Meta <- seurat_meta$meta[match(rownames(archr_meta), seurat_meta$atac_barcodes)]

proj <- addCellColData(proj, data = archr_meta$Isotype, name = "Isotype", cells = rownames(archr_meta))
proj <- addCellColData(proj, data = archr_meta$Subisotype, name = "Subisotype", cells = rownames(archr_meta))
proj <- addCellColData(proj, data = archr_meta$Meta, name = "Meta", cells = rownames(archr_meta))
saveArchRProject(proj)
