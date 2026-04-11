## Step 4: chromVAR-based motif enrichment and cross-method comparisons.
##
## Inputs:
## - MOCHA_SampleTileMatrix_Meta.rds
## - ChAI_PseudobulkRNA_Meta_Annotated.rds
## - Meta_EnrichedMotifs.csv (from Step 1)
##
## Outputs:
## - Meta_EnrichedMotifs_ChromVAR_ANOVA.csv
## - Meta_EnrichedMotifs_ChromVAR.csv
## - Meta_EnrichedMotifs_ChromVAR_Expressed.csv

library(MOCHA)
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(emmeans)

# Optionally load local development copy of MOCHA.
devtools::load_all("MOCHA")

# Load local helper implementations copied from ChAI.
helper_candidates <- c("helper_functions.r", "../helper_functions.r")
helper_path <- helper_candidates[file.exists(helper_candidates)][1]
if (is.na(helper_path)) {
  stop("Could not find helper_functions.r. Run from project root or analysis directory.")
}
source(helper_path)

run_dir <- "Run_Metatype_threshold"

STM <- readRDS(file.path(run_dir, "MOCHA_SampleTileMatrix_Meta.rds"))

chromSE <- makeChromVAR(STM, motifName = "Motifs")
saveRDS(chromSE, file.path(run_dir, "ChAI_ChromVAR_scores_Meta.rds"))

chromSE <- readRDS(file.path(run_dir, "ChAI_ChromVAR_scores_Meta.rds"))
chromSE <- chromSE$Z_Score
flatChrom <- flattenChAI(chromSE)

z_mat <- assay(flatChrom, "counts")
meta <- as.data.frame(colData(flatChrom))

# One-way cell-type ANOVA while controlling for donor/sample effect.
results <- lapply(seq_len(nrow(z_mat)), function(i) {
  df <- data.frame(
    z = z_mat[i, ],
    cell_type = meta$CellType,
    sample = sub("^.*__", "", meta$Sample)
  )
  fit <- lm(z ~ cell_type + sample, data = df)
  a <- anova(fit)

  data.frame(
    motif = rownames(z_mat)[i],
    p_value = a["cell_type", "Pr(>F)"]
  )
}) %>%
  bind_rows() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    mlogadjpval = -log10(pmax(p_adj, .Machine$double.xmin))
  )

get_ct_contrasts <- function(z, meta_df) {
  df <- data.frame(
    z = z,
    cell_type = meta_df$CellType,
    sample = sub("^.*__", "", meta_df$Sample)
  )
  fit <- lm(z ~ cell_type + sample, data = df)
  em <- emmeans(fit, ~ cell_type)
  as.data.frame(contrast(em, "eff"))
}

ct_contrasts <- lapply(results$motif, function(motif_name) {
  idx <- match(motif_name, rownames(z_mat))
  res <- get_ct_contrasts(z_mat[idx, ], meta)
  res$motif <- motif_name
  res
}) %>%
  bind_rows() %>%
  mutate(
    p_adj = p.adjust(p.value, "BH"),
    mlogadjpval = -log10(pmax(p_adj, .Machine$double.xmin))
  )

write.csv(ct_contrasts, file.path(run_dir, "Meta_EnrichedMotifs_ChromVAR_ANOVA.csv"), row.names = FALSE)

test_ct_vs_rest <- function(z_counts, meta_df, target_ct) {
  meta_df <- meta_df[colnames(z_counts), ]
  df <- meta_df %>% mutate(SampleID = sub("^.*__", "", Sample), cell_type = CellType)

  results_ct <- lapply(seq_len(nrow(z_counts)), function(i) {
    df$z <- z_counts[i, ]

    wide <- pivot_wider(df, id_cols = SampleID, names_from = cell_type, values_from = z)
    other_cols <- setdiff(colnames(wide), c("SampleID", target_ct))

    if (!(target_ct %in% colnames(wide)) || length(other_cols) == 0) {
      return(data.frame(motif = rownames(z_counts)[i], cell_type = target_ct, median_delta_z = NA, p_value = NA))
    }

    x <- wide[[target_ct]]
    y <- apply(wide[, other_cols], 1, mean, na.rm = TRUE)
    keep <- !is.na(x) & !is.na(y)

    if (sum(keep) < 2) {
      pval <- NA
      med_delta <- NA
    } else {
      pval <- wilcox.test(x[keep], y[keep], paired = TRUE, exact = FALSE)$p.value
      med_delta <- median(x[keep] - y[keep])
    }

    data.frame(motif = rownames(z_counts)[i], cell_type = target_ct, median_delta_z = med_delta, p_value = pval)
  })

  bind_rows(results_ct)
}

cell_types <- unique(meta$CellType)
results_all <- bind_rows(lapply(cell_types, function(ct) test_ct_vs_rest(z_mat, meta, ct))) %>%
  group_by(cell_type) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    mlogadjpval = -log10(pmax(p_adj, .Machine$double.xmin))
  ) %>%
  ungroup()

summary_signif <- results_all %>%
  filter(!is.na(p_adj), p_adj < 0.1, median_delta_z > 0) %>%
  group_by(cell_type) %>%
  summarise(n_significant_motifs = n(), .groups = "drop") %>%
  arrange(desc(n_significant_motifs))

summary_signif
write.csv(results_all, file.path(run_dir, "Meta_EnrichedMotifs_ChromVAR.csv"), row.names = FALSE)

wilcox_motifs <- results_all %>%
  filter(!is.na(p_adj), p_adj < 0.1, median_delta_z > 1)

# Compare chromVAR-enriched motifs with motif enrichment from Step 1.
chromModel <- model_ChromVAR(
  flatChrom,
  cellPopulation = "counts",
  modelFormula = exp ~ CellType + FragmentCounts + (1 | Sample),
  cellTypeMarkers = TRUE,
  numCores = 5
)
saveRDS(chromModel, file.path(run_dir, "MetaCluster_ChromModel.rds"))

chromModel <- readRDS(file.path(run_dir, "MetaCluster_ChromModel.rds"))
chromMotifs <- getEstimates(chromModel, "CellType_Marker")
table(chromMotifs$CellType)

allEnrichments <- read.csv(file.path(run_dir, "Meta_EnrichedMotifs.csv"))
dplyr::group_by(allEnrichments, CellType) %>%
  dplyr::filter(enrichment > 1) %>%
  dplyr::summarize(EnrichedMotifs = sum(adjp_val < 0.1))

motifListOriginal <- dplyr::group_by(allEnrichments, CellType) %>%
  dplyr::filter(enrichment > 1 & adjp_val < 0.1)

length(intersect(motifListOriginal$TFs, wilcox_motifs$motif))

overlap_per_celltype <- sapply(cell_types, function(ct) {
  sig_motifs <- wilcox_motifs %>%
    filter(cell_type == ct) %>%
    pull(motif)

  ref_motifs <- motifListOriginal %>%
    filter(CellType == ct) %>%
    pull(TFs)

  length(intersect(sig_motifs, ref_motifs))
})

overlap_per_celltype_df <- data.frame(
  cell_type = names(overlap_per_celltype),
  n_overlap = as.numeric(overlap_per_celltype)
)
overlap_per_celltype_df

table(motifListOriginal$CellType)

# Restrict motifs to TFs with expression support in matched cell type.
rnaSE <- readRDS(file.path(run_dir, "ChAI_PseudobulkRNA_Meta_Annotated.rds"))
exprGenes <- thresholdGenes(
  rnaSE,
  cellPopulation = "all",
  detectionThreshold = 0.001,
  expressionThreshold = 2,
  cellCountThreshold = 20
)
lengths(exprGenes)

motifList <- dplyr::group_by(results_all, cell_type) %>%
  dplyr::filter(median_delta_z > 0 & p_adj < 0.1)

for (ct in names(exprGenes)) {
  motifList <- dplyr::filter(
    motifList,
    (motif %in% exprGenes[[ct]] & cell_type == ct) | cell_type != ct
  )
}

table(motifList$cell_type)
write.csv(motifList, file.path(run_dir, "Meta_EnrichedMotifs_ChromVAR_Expressed.csv"), row.names = FALSE)
