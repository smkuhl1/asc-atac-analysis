## Step 1 (scRNA): Build pseudobulk RNA objects and annotate genes to hg38.
##
## Input:
## - teaseq_4_05.rds
##
## Outputs:
## - ChAI_PseudobulkRNA_Isotype.rds
## - ChAI_PseudobulkRNA_Isotype_Annotated.rds
## - ChAI_PseudobulkRNA_Meta.rds
## - ChAI_PseudobulkRNA_Meta_Annotated.rds

library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
library(Seurat)

# Load local helper implementations copied from ChAI.
helper_candidates <- c("helper_functions.r", "../helper_functions.r")
helper_path <- helper_candidates[file.exists(helper_candidates)][1]
if (is.na(helper_path)) {
    stop("Could not find helper_functions.r. Run from project root or analysis directory.")
}
source(helper_path)

seurat_rds <- "teaseq_4_05.rds"
seur <- readRDS(seurat_rds)

celltype_configs <- list(
    list(cellTypeColumn = "isotype", label = "Isotype"),
    list(cellTypeColumn = "meta", label = "Meta")
)

for (cfg in celltype_configs) {
    rnaSE <- makePseudobulkRNA(
        SO = seur,
        cellTypeColumn = cfg$cellTypeColumn,
        sampleColumn = "pbmc_sample_id",
        cellPopulations = "All",
        numCores = 1,
        Seurat_format = "SYMBOL",
        TxDb = NULL,
        OrgDb = NULL,
        normalize = TRUE
    )

    # Save unannotated object first because genome linking can drop weakly annotated genes.
    base_file <- paste0("ChAI_PseudobulkRNA_", cfg$label, ".rds")
    annot_file <- paste0("ChAI_PseudobulkRNA_", cfg$label, "_Annotated.rds")
    saveRDS(rnaSE, base_file)

    rnaSE_annot <- linkToGenome(
        rnaSE,
        TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
        OrgDb = org.Hs.eg.db
    )
    saveRDS(rnaSE_annot, annot_file)
}




