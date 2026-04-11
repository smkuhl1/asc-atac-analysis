## Annotate differential accessibility tiles with gene and region context.
##
## Inputs:
## - Run_Metatype_Threshold/MetaCluster_RegionModel.rds
##
## Output:
## - Run_Metatype_Threshold/DAT_table_annotated.csv

library(MOCHA)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# Optionally load local development copy of MOCHA.
devtools::load_all("MOCHA")

run_dir <- "Run_Metatype_Threshold"
region_model_rds <- file.path(run_dir, "MetaCluster_RegionModel.rds")
annotated_csv <- file.path(run_dir, "DAT_table_annotated.csv")

fix_ranges <- function(ranges, max_hg38_chr_size = 2.5e8) {
    parts <- do.call(rbind, strsplit(ranges, "[:-]"))
    chr <- parts[, 1]
    start <- as.numeric(parts[, 2])
    end_raw <- parts[, 3]

    bad <- which(as.numeric(end_raw) > max_hg38_chr_size)
    end_fixed <- end_raw
    if (length(bad) > 0) {
        end_fixed[bad] <- substr(end_raw[bad], 1, nchar(end_raw[bad]) - 1)
    }

    paste0(chr, ":", start, "-", end_fixed)
}

tileModel <- readRDS(region_model_rds)
continuousTiles <- getEstimates(tileModel, "CellType_Marker")
ziTiles <- getEstimates(tileModel, "ZI_CellType_Marker")
allTiles <- rbind(continuousTiles, ziTiles)
table(allTiles$CellType)

fixed_ranges <- fix_ranges(allTiles$Tiles)
gr <- StringsToGRanges(fixed_ranges)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- genes(txdb)
promoters_gr <- promoters(genes_gr, upstream = 2000, downstream = 200)
exons_gr <- exons(txdb)
introns_gr <- unlist(intronsByTranscript(txdb, use.names = TRUE))

ov_prom <- findOverlaps(gr, promoters_gr)
ov_gene <- findOverlaps(gr, genes_gr)
ov_exon <- findOverlaps(gr, exons_gr)
ov_intron <- findOverlaps(gr, introns_gr)

gr$promoter_gene_id <- NA_character_
gr$promoter_gene_id[queryHits(ov_prom)] <- as.character(promoters_gr$gene_id[subjectHits(ov_prom)])

gr$promoter_symbol <- mapIds(
    org.Hs.eg.db,
    keys = gr$promoter_gene_id,
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
)

gr$gene_id <- NA_character_
gr$gene_id[queryHits(ov_gene)] <- as.character(genes_gr$gene_id[subjectHits(ov_gene)])

gr$gene_symbol <- mapIds(
    org.Hs.eg.db,
    keys = gr$gene_id,
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
)

gr$annotation <- "intergenic"
gr$annotation[!is.na(gr$promoter_gene_id)] <- "promoter"
gr$annotation[queryHits(ov_exon)] <- "exon"
gr$annotation[queryHits(ov_intron)] <- "intron"

df <- cbind(allTiles, as.data.frame(mcols(gr)))

# Flatten list columns for CSV safety.
list_cols <- sapply(df, is.list)
if (any(list_cols)) {
    df[list_cols] <- lapply(df[list_cols], function(x) sapply(x, paste, collapse = ","))
}

write.csv(df, annotated_csv, row.names = FALSE)

