library(ChAI) ## from here: https://github.com/aifimmunology/ChAi
## devtools::load_all('ChAI') ## Load in the latest version. You can also use devtools::install('ChAI') to install it locally
library(MOCHA)
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db) 

devtools::load_all('MOCHA')

devtools::load_all('ChAI')
#### Analysis linking RNA and ATAC fo


setwd('Run_Metatype_threshold')
setwd('Run_Metatype_4v12')

STM <- readRDS('MOCHA_SampleTileMatrix_Meta.rds')

chromSE <- makeChromVAR(STM, motifName = 'Motifs')
saveRDS(chromSE, 'ChAI_ChromVAR_scores_Meta.rds')

chromSE <- readRDS('ChAI_ChromVAR_scores_Meta.rds')
chromSE <- chromSE$Z_Score
flatChrom <- flattenChAI(chromSE)

library(dplyr)
library(tidyr)

z <- assays(flatChrom)$counts   # motifs × (cell_type × sample)
meta <- as.data.frame(colData(flatChrom))



##########################################################
# ANOVA test
z_mat <- assay(flatChrom, "counts")
meta  <- as.data.frame(colData(flatChrom))

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


sig_motifs <- results %>%
  filter(p_adj < 0.05)

library(emmeans)

get_ct_contrasts <- function(z, meta) {
  df <- data.frame(
    z = z,
    cell_type = meta$CellType,
    sample = sub("^.*__", "", meta$Sample)
  )
  fit <- lm(z ~ cell_type + sample, data = df)
  em <- emmeans(fit, ~ cell_type)
  as.data.frame(contrast(em, "eff"))  # vs grand mean
}

ct_contrasts <- lapply(results$motif, function(m) {
  i <- match(m, rownames(z_mat))
  res <- get_ct_contrasts(z_mat[i, ], meta)
  res$motif <- m
  res
}) %>%
  bind_rows() %>%
  mutate(p_adj = p.adjust(p.value, "BH"),
        mlogadjpval = -log10(pmax(p_adj, .Machine$double.xmin)))

write.csv(ct_contrasts, 'Meta_EnrichedMotifs_ChromVAR_ANOVA.csv')










##################################
library(dplyr)
library(tidyr)

test_ct_vs_rest <- function(z_mat, meta, target_ct) {
  
  # Ensure metadata matches assay columns
  meta <- meta[colnames(z_mat), ]
  
  # Add motif values to metadata
  df <- meta %>%
    mutate(SampleID = sub("^.*__", "", meta$Sample), cell_type = CellType)
  
  # For each motif, compute median delta and paired Wilcoxon p-value
  results <- lapply(1:nrow(z_mat), function(i) {
    df$z <- z_mat[i, ]
    
    wide <- pivot_wider(df, id_cols = SampleID,
                        names_from = cell_type,
                        values_from = z)
    
    other_cols <- setdiff(colnames(wide), c("SampleID", target_ct))
    if (!(target_ct %in% colnames(wide)) || length(other_cols) == 0) {
      return(data.frame(motif = rownames(z_mat)[i],
                        cell_type = target_ct,
                        median_delta_z = NA,
                        p_value = NA))
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
    
    data.frame(motif = rownames(z_mat)[i],
               cell_type = target_ct,
               median_delta_z = med_delta,
               p_value = pval)
  })
  
  # Combine all motifs
  bind_rows(results)
}


################# run for all 5 cell types
cell_types <- unique(meta$CellType)

cell_types <- unique(meta$CellType)
results_all <- bind_rows(lapply(cell_types, function(ct) test_ct_vs_rest(z, meta, ct))) %>%
  group_by(cell_type) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"), 
        mlogadjpval = -log10(pmax(p_adj, .Machine$double.xmin))
        ) %>%
  ungroup()

summary_signif <- results_all %>%
  filter(!is.na(p_adj)) %>%       # drop any remaining NAs
  filter(p_adj < 0.1, median_delta_z > 0) %>%
  group_by(cell_type) %>%
  summarise(
    n_significant_motifs = n()
  ) %>%
  arrange(desc(n_significant_motifs))

summary_signif    

write.csv(results_all, 'Meta_EnrichedMotifs_ChromVAR.csv')
                                
wilcox_motifs <- results_all %>%
  filter(!is.na(p_adj)) %>%       # drop any remaining NAs
  filter(p_adj < 0.1, median_delta_z > 1)

####################################### linear modeling approach
chromModel = model_ChromVAR(flatChrom, 
            cellPopulation = 'counts',
            modelFormula = exp ~  CellType + FragmentCounts + (1|Sample),
            cellTypeMarkers = TRUE,
            numCores = 5)
saveRDS(chromModel, 'MetaCluster_ChromModel.rds')
                                
chromModel <- readRDS('MetaCluster_ChromModel.rds')
chromMotifs <- getEstimates(chromModel, 'CellType_Marker')
table(chromMotifs$CellType)

########################################
allEnrichments <- read.csv('Meta_EnrichedMotifs.csv')                                
dplyr::group_by(allEnrichments, CellType) %>% 
    dplyr::filter(enrichment > 1) %>%
    dplyr::summarize(EnrichedMotifs = sum(adjp_val < 0.1))

motifList = dplyr::group_by(allEnrichments, CellType) %>% 
    dplyr::filter(enrichment > 1 & adjp_val < 0.1) 

#######################################
length(intersect(motifList$TFs, wilcox_motifs$motif))                               

overlap_per_celltype <- sapply(cell_types, function(ct) {
  # Significant motifs for this cell type
  sig_motifs <- motifList %>%
    filter(cell_type == ct) %>%
    pull(motif)
  
  # Reference motifs for this cell type only
  ref_motifs <- motifListOriginal %>%
    filter(CellType == ct) %>%
    pull(TFs)
  
  # Count the overlap
  length(intersect(sig_motifs, ref_motifs))
})

# Convert to a data.frame
overlap_per_celltype_df <- data.frame(
  cell_type = names(overlap_per_celltype),
  n_overlap = as.numeric(overlap_per_celltype)
)

overlap_per_celltype_df


table(motifList$CellType)
table(motifListOriginal$CellType)


########################################

# filter expressed
rnaSE <- readRDS('/home/workspace/ChAI_PseudobulkRNA_Meta_Annotated.rds')
## A gene must be detected in at least 10% of cells for all samples with at least 10 cells to be counted as expressed
exprGenes = thresholdGenes(rnaSE,
                           cellPopulation = 'all', 
                           detectionThreshold = 0.001, 
                           expressionThreshold = 2,
                           cellCountThreshold = 20)
lengths(exprGenes)

## Let's get the significant motifs
motifList = dplyr::group_by(results_all, cell_type) %>% 
    dplyr::filter(median_delta_z > 0 & p_adj < 0.1) 

for(i in names(exprGenes)){

    motifList = dplyr::filter(motifList, 
                              (motif %in% exprGenes[[i]] & cell_type == i) |
                              cell_type != i)

}
table(motifList$cell_type) #Only TFs in ASC1 and ASC2?

write.csv(motifList, 'Meta_EnrichedMotifs_ChromVAR_Expressed.csv')







                                
                                