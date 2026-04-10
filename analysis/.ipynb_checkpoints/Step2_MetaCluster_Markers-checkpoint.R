library(ChAI) ## from here: https://github.com/aifimmunology/ChAi
## devtools::load_all('ChAI') ## Load in the latest version. You can also use devtools::install('ChAI') to install it locally
library(MOCHA)
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db) 

devtools::load_all('MOCHA')

devtools::load_all('ChAI')
#### Analysis linking RNA and ATAC for meta clusters

setwd('Run_Metatype_threshold')

system(sprintf("prlimit --pid=%d --as=50000000000", Sys.getpid()))



## Let's load in the data
STM <- readRDS('MOCHA_SampleTileMatrix_Meta.rds')
rnaSE <- readRDS('ChAI_PseudobulkRNA_Meta_Annotated.rds')

###
library(GSVAdata)
data(c2BroadSets)
geneSet1 = c2BroadSets[grep("^REACTOME", names(c2BroadSets))]


gseaSE <- runSSGSEA(rnaSE, geneSet = geneSet1, numCores= 1)
saveRDS(gseaSE, 'Meta_ssgsea.rds')

## Convert scATAC data to pseudobulk chromVAR data
chromSE <- makeChromVAR(STM, motifName = 'Motifs')
saveRDS(chromSE, 'ChAI_ChromVAR_scores_Meta.rds')
## pull out z-scores specifically
chromSE <- readRDS('ChAI_ChromVAR_scores_Meta.rds')
chromSE <- chromSE$Z_Score

## We then want to flatten it, so that rather than having one matrix for regions/gene by each plasma cell type,
## we instead have one matrix where each column is a combination of a cell type in one sample
## We will then model gene-tile relationships across celltypes and samples, controlling for sample-specific variation.

### Flat each data type
flatSTM <- combineSampleTileMatrix(STM)
flatRNA <- flattenChAI(rnaSE)
flatChrom <- flattenChAI(chromSE)
## Create metadata for aligning samples across modalities
flatRNA$SampleName = paste(flatRNA$CellType, flatRNA$pbmc_sample_id, sep= '_')
flatSTM$SampleName = paste(flatSTM$CellType, gsub(".*_","",flatSTM$Sample), sep= '_')
flatChrom$SampleName = paste(flatChrom$CellType, gsub(".*_","",flatChrom$Sample), sep= '_')

all(flatSTM$SampleName %in% flatRNA$SampleName) #TRUE


### get some DEGs so that we're not testing all genes near each tile, just ones that are significantly different
geneModel = model_scRNA(flatRNA, 
            cellPopulation = 'counts',
            modelFormula = exp ~  CellType + CellCounts + (1|pbmc_sample_id),
            ziFormula = ~ 0 + CellCounts,
            cellCountThreshold= 10,
            expressionThreshold = 2,
            cellTypeMarkers = TRUE,
            numCores = 25)
saveRDS(geneModel, 'MetaCluster_GeneModel.rds')

#Get DEGs
degList = getEstimates(geneModel, 'CellType_Marker')
table(degList$CellType)
write.csv(degList, 'MetaCluster_GeneModel_DEGs.csv')

### Repeat for scATAC-seq markers
STM <- readRDS('MOCHA_SampleTileMatrix_Meta.rds')
flatSTM <- combineSampleTileMatrix(STM)
regionModel = model_scATAC(flatSTM, 
            cellPopulation = 'counts',
            modelFormula = exp ~  CellType + FragmentCounts + (1|Sample),
            ziFormula = ~ FragmentCounts + CellType,
            cellTypeMarkers = TRUE,
            numCores = 40)
saveRDS(regionModel, 'MetaCluster_RegionModel.rds')

regionModel <- readRDS( 'MetaCluster_RegionModel.rds')

## let's look at these
continuousTiles = getEstimates(regionModel, 'CellType_Marker', FDR_threshold = 0.1)
## There are two portions in this test, because of zero-inflation
### 1. tests for the continuous changes: CellType_Markers
### 2. tests for changes in the number of empty values: ZI_CellType_Marker
### This is unique to scATAC, due to sparsity
ziTiles = getEstimates(regionModel, 'ZI_CellType_Marker', FDR_threshold = 0.1)
allTiles = rbind(continuousTiles, ziTiles)

table(allTiles$CellType) ## very few unique tiles for ASC3
## Let's get the matrix for ATAC accessibility
markerMat = pheatmap:::scale_rows(assays(flatSTM[unique(continuousTiles$Tiles),])[['counts']])
markerMat2 = pheatmap:::scale_rows(assays(flatSTM[unique(ziTiles$Tiles),])[['counts']])
markerMat3 =  pheatmap:::scale_rows(assays(flatSTM[unique(allTiles$Tiles),])[['counts']])

################################

markerMat = assays(flatSTM[unique(continuousTiles$Tiles),])[['counts']]
pdf('figures/test.pdf')
Heatmap(markerMat, name = 'Accessibility', top_annotation = column_annot, 
        cluster_rows= FALSE, cluster_columns=FALSE, 
        show_column_names = FALSE, show_row_names= FALSE, 
        use_raster = TRUE)
dev.off()

zero_mask <- markerMat == 0  # TRUE where original value was zero
scaled_mat <- pheatmap:::scale_rows(markerMat)
scaled_mat_with_marker <- scaled_mat
scaled_mat_with_marker[zero_mask] <- -1  # Must be outside the normal range [0, 1]
# Add a special color for -1 values
my_colors <- c("black", colorRampPalette(c("white", "blue", "red"))(100))  # black = zero peak

# Adjust breaks: include -1
breaksList <- c(-1, seq(0, 1, length.out = 101))

pdf('figures/test.pdf')

pheatmap(
  scaled_mat_with_marker,
  color = my_colors,
  breaks = breaksList,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Log2-scaled expression with zero peaks in black"
)





################################

### Extract metadata
atacMeta = as.data.frame(colData(flatSTM))

## Let's plot it now!
library(ComplexHeatmap)

metas <- paste0("ASC", 1:5)
meta.palette <- c("#1C7987", "#7286AC", "#7A3859", "#99C2A2", "#F0A7A0") %>% setNames(metas)

cluster_colors <- c(
  "ASC1" = "#1C7987",
  "ASC2" = "#7286AC",
  "ASC3" = "#7A3859",
  "ASC4" = "#99C2A2",
  "ASC5" = "#F0A7A0"
)

column_annot = HeatmapAnnotation(Cluster = atacMeta$CellType,
                                col = list(Cluster = cluster_colors)
                                 )




pdf('figures/MetaCluster_ChromatinHeatmap.pdf')
Heatmap(markerMat, name = 'Accessibility', top_annotation = column_annot, 
        cluster_rows= FALSE, cluster_columns=TRUE, 
        show_column_names = FALSE, show_row_names= FALSE, 
        use_raster = TRUE)
Heatmap(markerMat2, name = 'Accessibility', top_annotation = column_annot, 
        cluster_rows= FALSE, cluster_columns=TRUE, 
        show_column_names = FALSE, show_row_names= FALSE, 
        use_raster = TRUE)
                           
Heatmap(markerMat3, name = 'Accessibility', top_annotation = column_annot, 
        cluster_rows= FALSE, cluster_columns=TRUE, 
        show_column_names = FALSE, show_row_names= FALSE, 
        use_raster = TRUE)
dev.off()

### Repeat for chromVAR markers
chromModel = model_ChromVAR(flatChrom, 
            cellPopulation = 'counts',
            modelFormula = exp ~  CellType + FragmentCounts + (1|Sample),
            cellTypeMarkers = TRUE,
            numCores = 5)
saveRDS(chromModel, 'MetaCluster_ChromModel.rds')

## How many motifs?
## Let's try an alternative approach: chromVAR markers
## Let's pull them in
chromModel <- readRDS('MetaCluster_ChromModel.rds')
chromMotifs <- getEstimates(chromModel, 'CellType_Marker')
table(chromMotifs$CellType) ## even worse
## Not all show up here too. It's worse!

uniqueCellTypes = assayNames(STM)
allDEGs <- lapply(uniqueCellTypes, function(XX){
            flatSTM$Other = ifelse(flatSTM$CellType %in% XX, XX, 'Other')
            tmp = getDifferentialAccessibleTiles(flatSTM,
                cellPopulation = 'counts',
                groupColumn = 'Other',
                foreground = XX,
                numCores = 10,
                background = 'Other')
            tmp$CellType = XX
            tmp
    })
saveRDS(allDEGs, 'MetaCluster_DEGs_getDifferentialAccessibility.rds')
moreTiles =unlist(lapply(allDEGs, function(ZZ) plyranges::filter(ZZ, FDR < 0.1)$Tile))
# moreTiles = lapply(allDEGs, function(ZZ) plyranges::filter(ZZ, FDR < 0.1))
    
markerMat4 = assays(flatSTM[unique(moreTiles),])[['counts']]

## Let's plot it now!
library(ComplexHeatmap)
column_annot = HeatmapAnnotation(Cluster = atacMeta$CellType)

pdf('figures/MetaCluster_ChromatinHeatmap_getDifferentialAccessibiltiy.pdf')
Heatmap(log10(markerMat4+1), name = 'Accessibility', top_annotation = column_annot, 
        cluster_rows= FALSE, cluster_columns=TRUE, 
        show_column_names = FALSE, show_row_names= FALSE, 
        use_raster = TRUE)
dev.off()





ranges <- dat_table$X

parts <- do.call(rbind, strsplit(ranges, "[:-]"))
chr  <- parts[,1]
start <- as.numeric(parts[,2])
end   <- parts[,3]  # keep as character so we can trim

# Identify entries where the end coordinate is obviously wrong
bad <- which(as.numeric(end) > 2.5e8)   # >250 Mb = impossible for human genome

# Remove the final trailing digit for those bad ends
end_fixed <- end
end_fixed[bad] <- substr(end[bad], 1, nchar(end[bad]) - 1)

# Reassemble corrected ranges
ranges_fixed <- paste0(chr, ":", start, "-", end_fixed)

# Convert to GRanges
library(GenomicRanges)
gr <- StringsToGRanges(ranges_fixed)

                         
library(GenomicFeatures)
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb)
promoters <- promoters(genes(txdb), upstream = 2000, downstream = 200)
ov_prom <- findOverlaps(gr, promoters)

gr$promoter_gene_id <- NA
gr$promoter_gene_id[queryHits(ov_prom)] <- promoters$gene_id[subjectHits(ov_prom)]

library(org.Hs.eg.db)

gr$promoter_symbol <- mapIds(
    org.Hs.eg.db,
    keys = gr$promoter_gene_id,
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
)

genes <- genes(txdb)
ov_gene <- findOverlaps(gr, genes)

gr$gene_id <- NA
gr$gene_id[queryHits(ov_gene)] <- genes$gene_id[subjectHits(ov_gene)]
gr$gene_symbol <- mapIds(org.Hs.eg.db,
                         keys = gr$gene_id,
                         keytype = "ENTREZID",
                         column = "SYMBOL",
                         multiVals = "first")
                         
exons <- exons(txdb)
introns <- intronsByTranscript(txdb, use.names = TRUE)
introns <- unlist(introns)


gr$annotation <- "intergenic"

# promoter
gr$annotation[!is.na(gr$promoter_gene_id)] <- "promoter"

# exons
ov_exon <- findOverlaps(gr, exons)
gr$annotation[queryHits(ov_exon)] <- "exon"

# introns
ov_intron <- findOverlaps(gr, introns)
gr$annotation[queryHits(ov_intron)] <- "intron"

)


df_base <- as.data.frame(df)


                         