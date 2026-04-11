# Analysis Workflow Overview

This folder contains a four-step scATAC-seq + scRNA-seq analysis pipeline for integrative multi-omic characterization of antibody-secreting cell (ASC) metatypes.

## Prerequisites

Before running any scripts, ensure the following R packages and data are available:

**R Packages:**
- `ArchR` – single-cell ATAC peak calling and quality control
- `MOCHA` – multi-omic accessibility analysis (install from: https://github.com/aifimmunology/MOCHA)
- `ChAI` – chromatin and accessibility integration (install from: https://github.com/aifimmunology/ChAi)
- `SummarizedExperiment`, `TxDb.Hsapiens.UCSC.hg38.refGene`, `org.Hs.eg.db`, `BSgenome.Hsapiens.UCSC.hg38`
- `chromVARmotifs`, `ComplexHeatmap`, `dplyr`, `tidyr`, `emmeans`, `GSVAdata`

**Input Data:**
- Arrow files (.arrow)—placed in `input/` subdirectory
- Pre-processed Seurat object: `teaseq_4_05.rds` (contains cell metadata with isotype, metatype, and sample info)
- ChAI pseudobulk RNA object: `ChAI_PseudobulkRNA_Meta_Annotated.rds` (expected in `Run_Metatype_threshold/`)

## Script Execution Order

### **Step 0: Initialize ArchR Project** (`Step0_ArchR.R`)
- **Purpose:** Load Arrow files, create an ArchR project, and attach metadata from the Seurat object.
- **Outputs:**
  - `Plasma_ArchR/` – ArchR project directory with cell barcodes and metatype annotations
- **Runtime:** ~5–15 minutes (depends on number of cells)

### **Step 1: Call Peaks and Run Motif Analyses** (`Step1_Calling_Peaks_in_scATAC.R`)
- **Purpose:** Call open chromatin regions (tiles) per metatype, build a SampleTileMatrix, add insertion bias correction, motif annotation, and extract motif footprints.
- **Key Stages:**
  1. Call open tiles per cell type using MOCHA
  2. Generate SampleTileMatrix across samples (threshold = 0.3)
  3. Add Tn5 insertion bias correction
  4. Annotate motifs (chromVAR human PWMs, window width = 7)
  5. Identify cell type-specific unique regions
  6. Run motif enrichment (foreground = unique regions; background = all other regions)
  7. Extract top 50 enriched motifs per cell type; compute footprints (memory-intensive—processed one at a time)
  8. Export coverage tracks and insertion footprint files for visualization
- **Outputs (in `Run_Metatype_threshold/`):**
  - `MOCHA_SampleTileMatrix_Meta.rds`
  - `Meta_EnrichedMotifs.csv`
  - `Meta_Unique_Regions.rds`
  - `Meta_Files/` – accessibility exports for genome browser
  - `Meta_Files/Motif_Footprints/` – individual motif footprint objects per cell type
- **Runtime:** 2–4 hours (parallelized with 20–40 cores where available)

### **Step 2: Identify RNA/ATAC/chromVAR Markers** (`Step2_MetaCluster_Markers.R`)
- **Purpose:** Build statistical models linking metatype identity to accessibility, gene expression, and chromatin factor activity profiles.
- **Key Stages:**
  1. Compute Reactome pathway ssGSEA signatures from RNA
  2. Convert scATAC accessibility to chromVAR motif activity scores
  3. Flatten each modality to sample-level matrices (one row = one metatype–sample combination)
  4. Fit generalized linear mixed models on RNA, accessibility, and chromVAR with cell type as the fixed effect and sample as a random intercept
  5. Extract metatype markers (continuous and zero-inflated components for scATAC sparsity)
  6. Generate heatmaps of marker regions colored by metatype
- **Outputs (in `Run_Metatype_threshold/`):**
  - `Meta_ssgsea.rds` – Reactome pathway signatures
  - `MetaCluster_GeneModel.rds`, `MetaCluster_GeneModel_DEGs.csv`
  - `MetaCluster_RegionModel.rds`, `MetaCluster_ChromModel.rds`
  - `MetaCluster_DEGs_getDifferentialAccessibility.rds`
  - `figures/MetaCluster_ChromatinHeatmap.pdf`, `MetaCluster_ChromatinHeatmap_getDifferentialAccessibility.pdf`
- **Runtime:** 1–2 hours (depends on model fitting complexity; parallelized with 25–40 cores)

### **Step 4: chromVAR Enrichment and Cross-Method Comparison** (`Step4_ChromVAR.r`)
- **Purpose:** Validate and compare motif enrichment signals across statistical methods (ANOVA, Wilcoxon paired test, linear mixed model) and check concordance with expressed transcription factors.
- **Key Stages:**
  1. Compute cell type ANOVA on chromVAR scores while controlling for sample effect
  2. Extract per-cell type contrasts (vs. grand mean) using emmeans
  3. Run paired Wilcoxon test (each metatype vs. all others) for each motif across samples
  4. Fit linear mixed model on chromVAR scores with metatype and sample controls
  5. Compare results to motif enrichment from Step 1 (overrepresentation-based)
  6. Restrict final motif list to TFs with detectable expression in the corresponding metatype
- **Outputs (in `Run_Metatype_threshold/`):**
  - `Meta_EnrichedMotifs_ChromVAR_ANOVA.csv` – ANOVA contrasts
  - `Meta_EnrichedMotifs_ChromVAR.csv` – Wilcoxon test results
  - `Meta_EnrichedMotifs_ChromVAR_Expressed.csv` – expression-filtered motifs
- **Runtime:** 1–2 hours (depends on number of motifs × number of metatypes)

## Output Structure

```
Run_Metatype_threshold/
├── MOCHA_*.rds                          # SampleTileMatrix and open tiles
├── Meta_*.rds                           # Unique regions, ssGSEA, etc.
├── Meta_EnrichedMotifs*.csv             # Motif enrichment tables
├── MetaCluster_*.rds                    # Statistical model objects
├── MetaCluster_*_DEGs.csv               # Marker gene/region tables
├── figures/
│   ├── MetaCluster_ChromatinHeatmap.pdf
│   └── MetaCluster_ChromatinHeatmap_getDifferentialAccessibility.pdf
├── Meta_Files/
│   ├── SampleSpecificFiles/             # Coverage tracks per sample
│   └── Motif_Footprints/                # Per-motif footprint objects
└── Meta_MOCHA/                          # MOCHA working files
```

## Key Parameters

| Script | Parameter | Value | Note |
|--------|-----------|-------|------|
| Step 0 | HiSE data archive | 12 IDs | Deduplication to 6 unique samples |
| Step 1 | Tile threshold | 0.3 | Minimum fraction of samples with ≥1 read |
| Step 1 | Motif window width | 7 bp | chromVAR parameter |
| Step 1 | Footprint window | 500 bp | Half-window around motif peak |
| Step 2 | Expression threshold | 2 | pseudobulk count minimum |
| Step 2 | Cell count threshold | 10 | Minimum cells per sample to include |
| Step 4 | p-value threshold | 0.05 (ANOVA), 0.1 (Wilcoxon) | Adjusted; see outputs for details |

## Memory and Compute Requirements

- **Total runtime:** ~4–8 hours (all steps sequentially)
- **Peak memory:** 50+ GB (Step 1 motif footprint extraction; set memory limit or comment out if constrained)
- **Parallelization:** Scales well to 20–40 cores on supported systems; reduce `numCores` on smaller machines

## Troubleshooting

1. **Missing dependencies:** Run `devtools::install_github("aifimmunology/MOCHA")` and `devtools::install_github("aifimmunology/ChAi")` if packages are not found.
2. **Path issues:** Ensure all scripts are run from the workspace root directory; paths are relative and expect `input/`, `Plasma_ArchR/`, and `Run_Metatype_threshold/` subdirectories.
3. **Memory errors:** Reduce `numCores` or skip Step 1 motif footprint extraction; set `force = FALSE` in footprint export functions to skip already-computed files.
4. **Barcodes not matching:** Check that Seurat object and Arrow metadata use consistent barcode formats (pool_id, sample_id, cell barcode join strings).

## Citation

Analysis relies on:
- *ArchR*: Granja, J.M., Corces, M.R., et al. (Nat. Genet. 2021)
- *MOCHA*: Available at https://github.com/aifimmunology/MOCHA
- *ChAI*: Available at https://github.com/aifimmunology/ChAi
