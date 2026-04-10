library(ArchR)
library(MOCHA)
library(SummarizedExperiment)

#### These are transcript and genome references that need to be installed
library(TxDb.Hsapiens.UCSC.hg38.refGene) # transcript locations
library(org.Hs.eg.db) # gene names across database
library(BSgenome.Hsapiens.UCSC.hg38) # Genome assembly information for sequence info

devtools::load_all('MOCHA') ## Load in the latest version. You can also use devtools::install('MOCHA') to install it locally

proj <- loadArchRProject('Plasma_ArchR')
# Call open tiles for individual samples/celltypes
cellPopulations <- unique(proj$Meta)
cellPopulations <- cellPopulations[!is.na(cellPopulations)]

openTiles = callOpenTiles(proj, cellPopLabel = 'Meta',
                          cellPopulations = cellPopulations,
                          TxDb = 'TxDb.Hsapiens.UCSC.hg38.refGene',
                          OrgDb = 'org.Hs.eg.db', 
                          outDir = 'Meta_MOCHA',
                          numCores= 35)
saveRDS(openTiles, 'MOCHA_openTileCalls_Meta.rds')

setwd('Run_Metatype_threshold')
openTiles <- readRDS('/home/workspace/MOCHA_openTileCalls_Meta.rds')

# Genereate sample tile matrix across celltypes/samples
STM <- getSampleTileMatrix( tileResults = openTiles,                   
                              cellPopulations = "ALL",
                              groupColumn = NULL,
                              threshold = 0.3,
                              numCores = 20,
                              verbose = TRUE
                            )
saveRDS(STM, 'MOCHA_SampleTileMatrix_Meta.rds')                      

## Add Insertion bias. This allows us to correct for Tn5 insertion bias when looking at data. 
STM <- addInsertionBias(STM, numCores =10, verbose = TRUE)
saveRDS(STM, 'MOCHA_SampleTileMatrix_Meta.rds')  
gc()
## Add motif locations
library(chromVARmotifs)
data(human_pwms_v2)
STM <- addMotifSet(
   SampleTileObj = STM,
   motifPWMs = human_pwms_v2,
   returnSTM = TRUE, motifSetName = "Motifs", w = 7)
saveRDS(STM, 'MOCHA_SampleTileMatrix_Meta.rds') 

### Identify peaks called unique to each cell type
### start by by pulling out each cell type
cellTypeList <- getCellTypes(STM)
### Unique regions will be called as true for only one 'celltype'
uniqueTiles = rowRanges(STM)[which(rowSums(as.data.frame(mcols(rowRanges(STM)))[, cellTypeList]) == 1)]
saveRDS(uniqueTiles, 'Meta_Unique_Regions.rds')

### These are the files to use for motifEnrichment by cell type
### Let's create list of unique regions for each cell type
cellTypeTiles = lapply(cellTypeList, function(XX){
                    uniqueTiles[mcols(uniqueTiles)[[XX]]]        
            })
names(cellTypeTiles) = cellTypeList

### Let's run motif enrichments. 
motifPositions = STM@metadata$Motifs
allEnrichments = do.call('rbind', lapply(cellTypeList, function(XX){
        foregroundList = cellTypeTiles[[XX]]
        ## Get the background set of tiles - i.e. all the regions open across plasma cells
        backgroundList = plyranges::filter_by_non_overlaps(rowRanges(STM), foregroundList)
        tmpMat = MotifEnrichment(Group1 = foregroundList, Group2 = backgroundList, 
                        motifPosList = motifPositions)
        tmpMat$CellType = XX
        tmpMat
        }))
write.csv(allEnrichments, 'Meta_EnrichedMotifs.csv')

allEnrichments <- read.csv('Meta_EnrichedMotifs.csv')

## Too many enriched motifs. Let's select the top 50 per cell type and just validate those
dplyr::group_by(allEnrichments, CellType) %>% 
            dplyr::filter(enrichment > 1) %>% 
            dplyr::summarize(EnrichedMotifs = sum(adjp_val < 0.01))

allEnrichments$Motif_name <- rownames(allEnrichments)
## Top50 cell types
top50 = dplyr::group_by(allEnrichments, CellType) %>% 
        dplyr::filter(enrichment > 1 & adjp_val < 0.01) %>% dplyr::arrange(mlog10Padj) %>%
        dplyr::slice_head(n=50)

## Let's extract motif footprints and save them. This is a very memory intensive step, so we're just going to do it one motif at a time, and save it as a list of objects each time.
for(i in cellTypeList){
    filter_top50 <- dplyr::filter(top50, CellType == i)
    test_motifs <- filter_top50$Motif_name
    for(y in test_motifs){

        tmp = motifFootprint(SampleTileObj = STM,
                              motifName = 'Motifs',
                              specMotif = y,
                              regions = NULL,
                              cellPopulations = i,
                              windowSize = 500,
                              normTn5 = TRUE,
                              smoothTn5 = 10,
                              groupColumn = NULL,
                              subGroups = NULL,
                              sampleSpecific = FALSE,
                              numCores = 10,
                               force = FALSE,
                              verbose = FALSE) 
        gc()
        saveRDS(tmp, paste('Meta_Files/Motif_Footprints/CellType_', i, '_Motif_', y, '.rds'))
    }
}


#################################################################################################################

#### Exporting files for visualizations

################################################################################################################

STM <- readRDS('MOCHA_SampleTileMatrix_Meta.rds')

### Export accessibility files
exportCoverage(
  SampleTileObject = STM,
  dir = 'Meta_Files',
  # type = TRUE,
  cellPopulations = "ALL",
  groupColumn = NULL,
  subGroups = NULL,
  sampleSpecific = FALSE,
  saveFile = TRUE,
  numCores = 10,
  verbose = FALSE
)
gc()

#### Export simple-specific accessibility files
exportCoverage(
  SampleTileObject = STM,
  dir = 'Meta_Files/SampleSpecificFiles',
  # type = TRUE,
  cellPopulations = "ALL",
  groupColumn = NULL,
  subGroups = NULL,
  sampleSpecific = TRUE,
  saveFile = TRUE,
  numCores = 2,
  verbose = FALSE
)


## Export insertion files
exportLocalFootprints(SampleTileObj = STM,
                         cellPopulation ='All',
                         outDir = 'Meta_Files',
                         windowSize = 10,
                         groupColumn = NULL,
                         subGroups = NULL,
                         sampleSpecific = FALSE,
                         normTn5 = TRUE,
                         force=FALSE,
                         verbose=FALSE,
                         numCores = 6
                        )

exportLocalFootprints(SampleTileObj = STM,
                         cellPopulation ='All',
                         outDir = 'Meta_Files/SampleSpecificFiles',
                         windowSize = 10,
                         groupColumn = NULL,
                         subGroups = NULL,
                         sampleSpecific = TRUE,
                         normTn5 = TRUE,
                         force=TRUE,
                         verbose=FALSE,
                         numCores = 1
                        )
