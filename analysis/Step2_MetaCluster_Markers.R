## Step 2: Identify ATAC metatype markers and produce summary heatmaps.
##
## Inputs (from Step 1):
## - MOCHA_SampleTileMatrix_Meta.rds
##
## Outputs:
## - Marker model objects/tables and heatmap PDFs in figures/

library(MOCHA)
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
library(GSVAdata)
library(ComplexHeatmap)
library(dplyr)

# Optionally load local development copies.
devtools::load_all("MOCHA")

# Load local helper implementations copied from ChAI.
helper_candidates <- c("helper_functions.r", "../helper_functions.r")
helper_path <- helper_candidates[file.exists(helper_candidates)][1]
if (is.na(helper_path)) {
        stop("Could not find helper_functions.r. Run from project root or analysis directory.")
}
source(helper_path)

run_dir <- "Run_Metatype_threshold"
figure_dir <- file.path(run_dir, "figures")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# Keep Linux memory cap behavior where available, but do not fail on unsupported OS.
try(system(sprintf("prlimit --pid=%d --as=50000000000", Sys.getpid())), silent = TRUE)

# Load input object.
STM <- readRDS(file.path(run_dir, "MOCHA_SampleTileMatrix_Meta.rds"))

# Flatten modalities to sample-level matrices.
flatSTM <- combineSampleTileMatrix(STM)

flatSTM$SampleName <- paste(flatSTM$CellType, gsub(".*_", "", flatSTM$Sample), sep = "_")

# scATAC marker model.
regionModel <- model_scATAC(
        flatSTM,
        cellPopulation = "counts",
        modelFormula = exp ~ CellType + FragmentCounts + (1 | Sample),
        ziFormula = ~FragmentCounts + CellType,
        cellTypeMarkers = TRUE,
        numCores = 40
)
saveRDS(regionModel, file.path(run_dir, "MetaCluster_RegionModel.rds"))

regionModel <- readRDS(file.path(run_dir, "MetaCluster_RegionModel.rds"))
continuousTiles <- getEstimates(regionModel, "CellType_Marker", FDR_threshold = 0.1)
ziTiles <- getEstimates(regionModel, "ZI_CellType_Marker", FDR_threshold = 0.1)
allTiles <- rbind(continuousTiles, ziTiles)
table(allTiles$CellType)

markerMat <- pheatmap:::scale_rows(assays(flatSTM[unique(continuousTiles$Tiles), ])[["counts"]])
markerMat2 <- pheatmap:::scale_rows(assays(flatSTM[unique(ziTiles$Tiles), ])[["counts"]])
markerMat3 <- pheatmap:::scale_rows(assays(flatSTM[unique(allTiles$Tiles), ])[["counts"]])

# Accessibility heatmaps.
atacMeta <- as.data.frame(colData(flatSTM))
cluster_colors <- c(
        "ASC1" = "#1C7987",
        "ASC2" = "#7286AC",
        "ASC3" = "#7A3859",
        "ASC4" = "#99C2A2",
        "ASC5" = "#F0A7A0"
)
column_annot <- HeatmapAnnotation(Cluster = atacMeta$CellType, col = list(Cluster = cluster_colors))

pdf(file.path(figure_dir, "MetaCluster_ChromatinHeatmap.pdf"))
Heatmap(
        markerMat,
        name = "Accessibility",
        top_annotation = column_annot,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        use_raster = TRUE
)
Heatmap(
        markerMat2,
        name = "Accessibility",
        top_annotation = column_annot,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        use_raster = TRUE
)
Heatmap(
        markerMat3,
        name = "Accessibility",
        top_annotation = column_annot,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        use_raster = TRUE
)
dev.off()


# Alternative differential accessibility approach.
uniqueCellTypes <- assayNames(STM)
allDEGs <- lapply(uniqueCellTypes, function(ct) {
        flatSTM$Other <- ifelse(flatSTM$CellType %in% ct, ct, "Other")
        tmp <- getDifferentialAccessibleTiles(
                flatSTM,
                cellPopulation = "counts",
                groupColumn = "Other",
                foreground = ct,
                numCores = 10,
                background = "Other"
        )
        tmp$CellType <- ct
        tmp
})
saveRDS(allDEGs, file.path(run_dir, "MetaCluster_DEGs_getDifferentialAccessibility.rds"))

moreTiles <- unlist(lapply(allDEGs, function(x) plyranges::filter(x, FDR < 0.1)$Tile))
markerMat4 <- assays(flatSTM[unique(moreTiles), ])[["counts"]]
column_annot <- HeatmapAnnotation(Cluster = atacMeta$CellType)

pdf(file.path(figure_dir, "MetaCluster_ChromatinHeatmap_getDifferentialAccessibility.pdf"))
Heatmap(
        log10(markerMat4 + 1),
        name = "Accessibility",
        top_annotation = column_annot,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        use_raster = TRUE
)
dev.off()


                         