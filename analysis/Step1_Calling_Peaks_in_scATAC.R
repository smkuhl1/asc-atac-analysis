## Step 1: Call open chromatin peaks by metatype and run motif-centric analyses.
##
## Inputs:
## - ArchR project: Plasma_ArchR/
##
## Outputs:
## - MOCHA open tiles, SampleTileMatrix, motif enrichment tables, motif footprints,
##   and coverage/footprint exports under the run directory.

library(ArchR)
library(MOCHA)
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)

# Optionally load local development copy of MOCHA.
devtools::load_all("MOCHA")

project_dir <- "Plasma_ArchR"
run_dir <- "Run_Metatype_threshold"
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

open_tiles_rds <- file.path(run_dir, "MOCHA_openTileCalls_Meta.rds")
stm_rds <- file.path(run_dir, "MOCHA_SampleTileMatrix_Meta.rds")
unique_regions_rds <- file.path(run_dir, "Meta_Unique_Regions.rds")
motif_enrichment_csv <- file.path(run_dir, "Meta_EnrichedMotifs.csv")
meta_files_dir <- file.path(run_dir, "Meta_Files")
footprint_dir <- file.path(meta_files_dir, "Motif_Footprints")
sample_specific_dir <- file.path(meta_files_dir, "SampleSpecificFiles")
dir.create(footprint_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(sample_specific_dir, showWarnings = FALSE, recursive = TRUE)

proj <- loadArchRProject(project_dir)
cell_populations <- unique(proj$Meta)
cell_populations <- cell_populations[!is.na(cell_populations)]

# Call open tiles per metatype.
openTiles <- callOpenTiles(
  proj,
  cellPopLabel = "Meta",
  cellPopulations = cell_populations,
  TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene",
  OrgDb = "org.Hs.eg.db",
  outDir = file.path(run_dir, "Meta_MOCHA"),
  numCores = 35
)
saveRDS(openTiles, open_tiles_rds)

openTiles <- readRDS(open_tiles_rds)

# Generate sample tile matrix across cell types and samples.
STM <- getSampleTileMatrix(
  tileResults = openTiles,
  cellPopulations = "ALL",
  groupColumn = NULL,
  threshold = 0.3,
  numCores = 20,
  verbose = TRUE
)
saveRDS(STM, stm_rds)

# Add insertion bias and motifs.
STM <- addInsertionBias(STM, numCores = 10, verbose = TRUE)
saveRDS(STM, stm_rds)
gc()

data(human_pwms_v2)
STM <- addMotifSet(
  SampleTileObj = STM,
  motifPWMs = human_pwms_v2,
  returnSTM = TRUE,
  motifSetName = "Motifs",
  w = 7
)
saveRDS(STM, stm_rds)

# Identify cell type-specific unique regions.
cell_type_list <- getCellTypes(STM)
unique_tiles <- rowRanges(STM)[
  which(rowSums(as.data.frame(mcols(rowRanges(STM)))[, cell_type_list]) == 1)
]
saveRDS(unique_tiles, unique_regions_rds)


STM <- readRDS(stm_rds)

# Export coverage tracks.
exportCoverage(
  SampleTileObject = STM,
  dir = meta_files_dir,
  cellPopulations = "ALL",
  groupColumn = NULL,
  subGroups = NULL,
  sampleSpecific = FALSE,
  saveFile = TRUE,
  numCores = 10,
  verbose = FALSE
)
gc()

exportCoverage(
  SampleTileObject = STM,
  dir = sample_specific_dir,
  cellPopulations = "ALL",
  groupColumn = NULL,
  subGroups = NULL,
  sampleSpecific = TRUE,
  saveFile = TRUE,
  numCores = 2,
  verbose = FALSE
)

# Export local insertion footprint tracks.
exportLocalFootprints(
  SampleTileObj = STM,
  cellPopulation = "All",
  outDir = meta_files_dir,
  windowSize = 10,
  groupColumn = NULL,
  subGroups = NULL,
  sampleSpecific = FALSE,
  normTn5 = TRUE,
  force = FALSE,
  verbose = FALSE,
  numCores = 6
)

exportLocalFootprints(
  SampleTileObj = STM,
  cellPopulation = "All",
  outDir = sample_specific_dir,
  windowSize = 10,
  groupColumn = NULL,
  subGroups = NULL,
  sampleSpecific = TRUE,
  normTn5 = TRUE,
  force = TRUE,
  verbose = FALSE,
  numCores = 1
)
