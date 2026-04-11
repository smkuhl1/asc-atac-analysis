# take select motifs that have differences across subsets and finalize for manuscript inclusion

library(MOCHA)
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db) 
library(doParallel)

devtools::load_all('/home/workspace/MOCHA')

# Load local helper implementations copied from ChAI.
helper_candidates <- c('helper_functions.r', '../helper_functions.r')
helper_path <- helper_candidates[file.exists(helper_candidates)][1]
if (is.na(helper_path)) {
  stop('Could not find helper_functions.r. Run from project root or figures directory.')
}
source(helper_path)

footprintList <- c('STAT2', 'IRF9', 'IRF3', 'IRF2', 'IRF1', 'NFIX')

footprintList <- c('NFIA')
footprintList <- c('BATF')
footprintList <- c('SP3', 'SP4', 'RUNX3')
footprintList <- c('RUNX3')


setwd('/home/workspace/Run_Metatype_threshold')
STM <- readRDS('MOCHA_SampleTileMatrix_Meta.rds')

# dat_tiles <- dat %>% filter(CellType == 'ASC2')
# calculate motifFootprint for each and save as RDS file
for(y in footprintList){
    print(y)
    tmp = motifFootprint(SampleTileObj = STM,
                          motifName = 'Motifs',
                          specMotif = c(y),
                          regions = NULL,
                          cellPopulations = c('ASC1', 'ASC2', 'ASC3', 'ASC4', 'ASC5'),
                          windowSize = 200,
                          normTn5 = TRUE,
                          smoothTn5 = 20,
                          groupColumn = NULL,
                          subGroups = NULL,
                          sampleSpecific = FALSE,
                          numCores = 15,
                           force = TRUE,
                          verbose = FALSE) 
    gc()
    saveRDS(tmp, paste0('/home/workspace/Motif_Footprints/footprint_plots_final/Meta_Footprint_', y, '_dat.rds'))
    rm(tmp)
}


# plot footprints and save pdf
for (y in footprintList) {
    print(y)
    imagePath <- paste0('/home/workspace/Motif_Footprints/footprint_plots_final/Plot_Meta_Footprint_', y, '_dat.pdf')
    footprintPath <- paste0('/home/workspace/Motif_Footprints/footprint_plots_final/Meta_Footprint_', y, '_dat.rds')
    
    # footprintPath <- '/home/workspace/Motif_Footprints/round2/Meta_Footprint_NFIA.rds'
    footprint <- readRDS(footprintPath)
    footprintName1 <- paste0('ASC1__', y)
    footprintName2 <- paste0('ASC2__', y)
    footprintName3 <- paste0('ASC3__', y)
    footprintName4 <- paste0('ASC4__', y)
    footprintName5 <- paste0('ASC5__', y)
        
    p <- plotMotifs(motifSE = footprint, 
               footprint = c(footprintName1, 
                             footprintName2, 
                             footprintName3,
                             footprintName4,
                             footprintName5
                            ), 
               groupColumn = NULL, 
               returnDF = FALSE, 
               returnPlotList = TRUE,
               topPercentage = 0.1,
               plotIndividualRegions = TRUE,
               relHeights = c(0.3, 0.7)
                              )
    print(unique(p$MotifAverage$data$PlotGroup))
    
    dt <- p$MotifAverage$data
    dt[, PlotGroup := sub(".*(ASC[1-5]).*", "\\1", PlotGroup)]
    p$MotifAverage$data <- dt
    print(unique(p$MotifAverage$data$PlotGroup))
    
    # change colors to match subsets
    p$MotifAverage <- p$MotifAverage +
      ggplot2::scale_color_manual(
        name = "ASC Subset",
        values = c(
        "ASC1" = "#1C7987",
        "ASC2" = "#7286AC",
        "ASC3" = "#7A3859", 
        "ASC4" = "#99C2A2",
        "ASC5" = "#F0A7A0"
      )) +
      ggplot2::theme(legend.position = 'right') +
      ggplot2::ggtitle(paste(y, "Motif Footprint"))
      
    
    ggplot2::ggsave(
      filename = imagePath,
      plot = p$MotifAverage,
      width = 6, height = 4, units = "in"
    )

}
















    
plots$MotifAverage

dev.off()

# recreate the combined figure
combined <- cowplot::plot_grid(
  plotlist = list(
    p$MotifAverage + ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ),
    p$LocationSpecific
  ),
  ncol = 1, align = "v", axis = "brl",
  rel_heights = c(0.3, 0.7)  # match what you want
)

grDevices::pdf(imagePath, width = 8, height = 10)
print(combined)
grDevices::dev.off()



