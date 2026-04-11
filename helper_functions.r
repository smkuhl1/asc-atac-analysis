# helper functions originally from ChAI 


#' @title High-throughut Modeling of scATAC data
#'
#' @description \code{model_scATAC} High-throughput modeling of tile accessibility from scATAC \code{\link[glmmTMB]{glmmTMB}}. 
#'
#' @param atacSE A SummarizedExperiment object generated from
#'   getSampleTileMatrix. 
#' @param cellPopulation Name of a cell type. 
#' @param modelFormula The formula for the continuous data that should be used within glmmTMB. It should be in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the atacSE metadata, except for CellType, FragmentCounts and CellCount, which will be extracted from the atacSE.
#'   modelFormula must start with 'exp' as the response.
#'   See \link[glmmTMB]{glmmTMB}.
#' @param ziFormula The formula for the zero-inflated data that should be used within glmmTMB. It should be in the
#'   format ( ~ factors). All factors must be found in column names
#'   of the atacSE colData metadata, except for CellType, FragmentCounts and CellCount, which will be extracted from the atacSE.
#' @param zi_threshold Zero-inflated threshold ( range = 0-1), representing the fraction of samples with zeros. At or above this threshold, the zero-inflated modeling kicks in.
#' @param initialSampling Size of data to use for pilot
#' @param cellTypeMarkers A boolean, used only if you want to identify cell type markers. 
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return results a SummarizedExperiment containing modeling results. Assays are metrics related to the model coefficients,
#'          including the Estimate, Std_Error, df, t_value, p_value. Within each assay, each row corresponds to each row of
#'          the SummarizedExperiment and columns correspond to each fixed effect variable within the model.
#'          Any row metadata from the ExperimentObject (see rowData(ExperimentObj)) is preserved in the output. 
#'          The Residual matrix and the variance of the random effects are saved in the metadata slot of the output.
#'          If multiple cell types are provided, then a named list of the above objects will be returned, one for each object. 
#'
#'
#'
#' @examples
#' \dontrun{
#'   modelList <- model_scATAC(STM[c(1:1000),], 
#'                  cellPopulation = 'CD16 Mono',
#'                  modelFormula = exp~ Age + Sex + days_since_symptoms + (1|PTID), 
#'                  ziFormula = ~ 0 + FragmentCounts + Age, 
#'                  verbose = TRUE, 
#'                  numCores = 35 )
#' }
#'
#' @export
#' @keywords modeling_individual
model_scATAC <- function(atacSE,
                      cellPopulation,
                      modelFormula = NULL,
                      ziFormula = ~ 0,
                      zi_threshold = 0,
                      initialSampling = 5,
                      cellTypeMarkers = FALSE,
                      verbose = FALSE,
                      numCores = 2) {
  if (!requireNamespace("MOCHA", quietly = TRUE)) {
      stop(
      "Package 'MOCHA' is required for pilot_scATAC. ",
      "Please install 'MOCHA' to proceed."
      )
  }

  if(isChAIObject(atacSE, type = 'data', returnType = TRUE) != 'scATAC'){

    stop('atacSE is not a MOCHA SampleTileObject.')

   }
    
  if (!methods::is(modelFormula,'formula')) {
    stop("modelFormula was not provided as a formula.")
  }
    
    if (!methods::is(ziFormula,'formula')) {
    stop("ziFormula was not provided as a formula.")
  }
      
  if (zi_threshold < 0 | zi_threshold > 1 | ! is.numeric(zi_threshold)) {
    stop("zi_threshold must be between 0 and 1.")
  }
    
    
  ## Evaluate cellPopulations. 
  if(all(tolower(cellPopulation) == 'all')){
   
      cellPopulation = SummarizedExperiment::assayNames(atacSE)
      
  }
    
      
  if(!all(cellPopulation %in% SummarizedExperiment::assayNames(atacSE))){
   
      stop('cellPopulation not found within atacSE object')
      
  }
      
  #If multiple cell types, then loop over all and return a list of model objects. 
  if(length(cellPopulation) > 1){
    
    allRes = lapply(cellPopulation, function(XX){
            message(stringr::str_interp('Modeling ${XX}'))
            tryCatch({
                gc()
                model_scATAC(atacSE,
                      cellPopulation = XX,
                      modelFormula = modelFormula,
                      ziFormula = ziFormula,
                      zi_threshold = zi_threshold,
                      initialSampling = initialSampling,
                      cellTypeMarkers = cellTypeMarkers,
                      verbose = verbose,
                      numCores = numCores)
                }, error = function(e){e})
        })
    names(allRes) = cellPopulation
    return(allRes)
    
  } else if(cellPopulation == 'counts' & !cellTypeMarkers){
      
    newObj <- atacSE

 }else if(cellPopulation == 'counts' & cellTypeMarkers){

     variableList = all.vars(modelFormula)
     if(!any(variableList %in% 'CellType')){
        stop('You have set cellTypeMarkers = TRUE. Please provided a flattened MOCHA object using combineSampleTileMatrix and include CellType as a variable in your formula.')
      }
     allCellTypes = unique(atacSE$CellType)
     mcolData = GenomicRanges::mcols(SummarizedExperiment::rowRanges(atacSE))

     if(!all(allCellTypes %in% colnames(mcolData))){
        stop('Misaligned cell type names. The CellType column contains names of cell types that MOCHA does not appear to have been called on. Please confirm that all cell types names are aligned between CellType column within the metadata, and the column names of the rowRanges slot of the MOCHA object. Each column name within the rowRanges slot records whether a region was called as open or not in a given cell type.')
     }
     LHS <- deparse(modelFormula)
     modelFormula = as.formula(gsub("CellType", "CellType_", LHS))

     ziF <- deparse(ziFormula)
     ziFormula = as.formula(gsub("CellType", "CellType_", ziF))
      
     allRes = lapply(allCellTypes, function(XX){
            message(stringr::str_interp('Identifying cell type markers for ${XX}'))
            tryCatch({
                gc()
                atacSE$CellType_ = ifelse(atacSE$CellType == XX, 'Marker', 'Background')
                atacSE$CellType_ = factor(atacSE$CellType_, levels = c('Background','Marker'))
                regList = mcolData[[XX]]
                model_scATAC(atacSE[regList,],
                      cellPopulation = 'counts', 
                        modelFormula = modelFormula, 
                        ziFormula = ziFormula,
                        zi_threshold = zi_threshold,
                      initialSampling = initialSampling,
                      cellTypeMarkers = FALSE,
                      verbose = verbose,
                      numCores = numCores)
                
            }, error = function(e){e})
        })
    names(allRes) = allCellTypes
    return(allRes)

  }else{
      
    newObj <- MOCHA::combineSampleTileMatrix(MOCHA::subsetMOCHAObject(atacSE, subsetBy = 'celltype', groupList = cellPopulation, subsetPeaks = TRUE))
      
  }

  exp <- .model_generic(SE_Object = newObj,
                      continuousFormula = modelFormula,
                      ziFormula = ziFormula,
                      zi_threshold = zi_threshold,
                      initialSampling = initialSampling,
                      family = stats::gaussian(),
                      modality = 'scATAC_Model',
                      verbose = verbose,
                      numCores = numCores)
  return(exp)

}




#' @title High-through modeling (Generalized)
#'
#' @description \code{model_General} Runs generalized GLM modeling for
#'   continuous, non-zero inflated data using \code{\link[glmmTMB]{glmmTMB}}
#'
#' @param ExperimentObj A SummarizedExperiment object generated from chromVAR, or other.
#'   It is expected to contain normally distributed data without zero-inflation. 
#' @param assayName The name of the assay to model within the SummarizedExperiment. 
#' @param modelFormula The formula to use, in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the ExperimentObj metadata. modelFormula must start with 'exp' as the response.
#'   See \link[glmmTMB]{glmmTMB}.
#' @param ziFormula A formula for zero inflation. By default, zero-inflated modeling is turned off by setting ziFormula = ~ 0
#' @param family distribution family parameter, passed to glmmTMB to describe the data's distribution.
#'     Default is normal (gaussian()). See  \link[glmmTMB]{glmmTMB}.
#' @param initialSampling Size of data to use for pilot
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return results a SummarizedExperiment containing modeling results. Assays are metrics related to the model coefficients,
#'          including the Estimate, Std_Error, df, t_value, p_value. Within each assay, each row corresponds to each row of
#'          the SummarizedExperiment and columns correspond to each fixed effect variable within the model.
#'          Any row metadata from the ExperimentObject (see rowData(ExperimentObj)) is preserved in the output. 
#'          The Residual matrix and the variance of the random effects are saved in the metadata slot of the output.
#'          If multiple cell types are provided, then a named list of the above objects will be returned, one for each object. 
#'
#'
#'
#' @examples
#' \dontrun{
#'   modelList <- model_General(ExperimentObj,
#'    assayName = "z",
#'     modelFormula = NULL,
#'     initialSampling = 5,
#'     verbose = FALSE,
#'     numCores = 1
#'  )
#' }
#'
#' @export
#' @keywords modeling_individual
model_General <- function(ExperimentObj,
                    assayName = 'General',
                    modelFormula,
                    ziFormula = ~ 0,
                    family = stats::gaussian(),
                    initialSampling = 5,
                    verbose = FALSE,
                    numCores = 2) {

  if(isChAIObject(ExperimentObj, type = 'data', returnType = TRUE) != 'General'){

    stop('ExperimentObj is not a ChAI General Modality object, generated via importGeneralModality.')

  }
    
   if (methods::is(ziFormula, 'character')) {
    ziFormula <- stats::as.formula(ziFormula)
  }

  
  if(!any(names(SummarizedExperiment::assays(ExperimentObj)) %in% assayName)){
    stop('ExperimentObj does not contain an assay that matches the assayName input variable.')
  }
  SummarizedExperiment::assays(ExperimentObj) = SummarizedExperiment::assays(ExperimentObj)[assayName]

  exp <- .model_generic(SE_Object = ExperimentObj,
                      continuousFormula = stats::as.formula(modelFormula),
                      ziFormula = ziFormula,
                      zi_threshold = 0,
                      family = family,
                      initialSampling = initialSampling,
                      modality = 'General_Model',
                      verbose = verbose,
                      numCores = numCores)
  return(exp)

}

#' @title Generate a ChAI ChromVAR Object from MOCHA's Tile-by-Sample Matrix Object
#'
#' @description \code{makeChromVAR} Runs ChromVAR on MOCHA's Tile-by-Sample Matrix object. Slow funciton. 
#'
#' @param atacSE A MOCHA Tile-by-Sample Object (SummarizedExperiment)generated from getSampleTileMatrix within \code{\link[MOCHA]{MOCHA}}. 
#' @param cellPopulation Names of cell types to analyze. Must match assay names in the SummarizedExperiment object. 
#'  Alternative, if you want to run ChromVAR across all celltypes, you can provide the output of combineSampleTileMatrix, and set this parameter to 'counts'.
#' @param motifName Name of metadata slot that has motif information for analysis.
#' @param withinCellType Boolean. Default is FALSE, in which case accessibility across celltypes will be merged into one matrix for ChromVAR, rather than running ChromVAR seperately on the matrix for each cell type when withinCellType = TRUE.
#' @param exportRaw Boolean. Default is FALSE, and will export a ChAI object with ChromVAR deviations reformated for modeling. If set to true, it will either export the original chromVARDeviations object for analysis, either as a list by cell type if withinCellType is TRUE or as a single object if withinCellType is TRUE. This chromVARDeviations object can be later reformatted for ChAI using reformatChromVAR.
#' @param numCores Default is 1. Uses chromVAR's standard parallelization, which has memory leak errors at times. Use at your own risk. 
#' @param verbose Boolean.
#' @return A named list of ChromVAR objects. 
#' 
#' @keywords data_import
#'
#' @export

makeChromVAR <- function(atacSE, motifName ='Motifs',
                      cellPopulation = NULL,
                      withinCellType = FALSE,
                      exportRaw = FALSE,
                    numCores = 1,
                      verbose = TRUE) {

    if (!requireNamespace("chromVAR", quietly = TRUE)) {
        
        stop(
        "Package 'chromVAR' is required for makeChromVAR. ",
        "Please install 'chromVAR' to proceed."
        )
        
    }
    
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
        
        stop(
        "Package 'BiocParallel' is required for makeChromVAR. ",
        "Please install 'BiocParallel' to proceed."
        )
        
    }

    if(is.null(cellPopulation)){

        cellPopulation = names(SummarizedExperiment::assays(atacSE))

    }

    if(is.null(cellPopulation)){

        cellPopulation = names(SummarizedExperiment::assays(atacSE))

    }

    if (
        any(!cellPopulation %in% names(SummarizedExperiment::assays(atacSE)))
    ) {

        stop("cellPopulation was not found within atacSE.")

    } else if (withinCellType) {
   
        #Generate a list of combined objects
        newObj <- lapply(cellPopulation, function(x){
             MOCHA::combineSampleTileMatrix(MOCHA::subsetMOCHAObject(atacSE, subsetBy = 'celltype', groupList = x, subsetPeaks = TRUE))
        })
        names(newObj) <- cellPopulation

    } else if(all(cellPopulation == 'counts')){
        newObj <- atacSE
    } else {
        newObj <- MOCHA::combineSampleTileMatrix(MOCHA::subsetMOCHAObject(atacSE, 
                                                                          subsetBy = 'celltype', groupList = cellPopulation, subsetPeaks = TRUE))
    }
    
    genome <- S4Vectors::metadata(atacSE)$Genome
    genome <- BSgenome::getBSgenome(genome)

    if(numCores > 1){
        
        BiocParallel::register(BiocParallel::SerialParam())
        
    }else if(Sys.info()[['sysname']] != 'Windows'){
        
        BiocParallel::register(BiocParallel::MulticoreParam(numCores, progressbar = TRUE))
        
    }else{
        
        BiocParallel::register(BiocParallel::SnowParam(workers = numCores, type = "SOCK"))
        
    }

    
    #Either iterate ovewr a list of newObj, or just directly on newObj to generate ChromVAR
    if(any(tolower(class(newObj)) %in% 'list')){

        chromVAROut <- lapply(cellPopulation, function(XX){
            if(verbose){ message('Analyzing ', XX)}
            newObj[[XX]] <- addGCBias_ChAI(newObj[[XX]], genome, verbose)
            anno_ix <- chromVAR::getAnnotations(newObj[[XX]]@metadata[[motifName]], 
                        rowRanges = SummarizedExperiment::rowRanges(newObj[[XX]]))
            chromVAR::computeDeviations(object =newObj[[XX]], 
                annotations = anno_ix)
        })
         names(chromVAROut) <- cellPopulation
        
        if(exportRaw){
            return(chromVAROut)
        }

    }else{

        if(verbose){ message('Analyzing ', cellPopulation)}
        newObj <- addGCBias_ChAI(newObj, genome, verbose)

        anno_ix <- chromVAR::getAnnotations(newObj@metadata[[motifName]], 
                        rowRanges = SummarizedExperiment::rowRanges(newObj))
        chromVAROut <- chromVAR::computeDeviations(object = newObj, 
                    annotations = anno_ix)
        
        if(exportRaw){
            return(chromVAROut)
        }
    
    }

    newOut_Dev <- reformatChromVAR(chromVAROut, selectDev =TRUE)
    newOut_Z <- reformatChromVAR(chromVAROut, selectDev =FALSE)

    newOut <- list('Z_Score' = newOut_Z, 'Deviations' = newOut_Dev)

    BiocParallel::register(BiocParallel::SerialParam())


    return(newOut)

}




#' @title \code{reformatChromVARList}
#'
#' @description \code{reformatChromVARList} pseodubulks a Seurat object by sample and cell type into a SummarizedExperiment, similar to MOCHA
#'
#' @param atacSE A SummarizedExperiment object generated from
#'   getSampleTileMatrix or combineSampleTileMatrix
#' @param genome BS Genome file to use as reference for GC bias scoring. 
#' @param verbose Boolean. verbose flag. Default is FALSE. 
#' @return A chromVAR object with GC bias
#'
#' @noRd

addGCBias_ChAI <- function(obj1, genome, verbose = FALSE){

    if (!requireNamespace("chromVAR", quietly = TRUE)) {
        stop(
        "Package 'chromVAR' is required for addGCBias_ChAI. ",
        "Please install 'chromVAR' to proceed."
        )
    }

    obj1 <- chromVAR::addGCBias(obj1, genome = genome)
    if (any(is.na(SummarizedExperiment::rowData(obj1)$bias))) {
        naList <- is.na(SummarizedExperiment::rowData(obj1)$bias)
        
        if (verbose) {
        warning(paste(sum(naList), "NaNs found within GC Bias", sep = " "))
        }
        
        SummarizedExperiment::rowData(obj1)$bias[which(naList)] <- mean(SummarizedExperiment::rowData(obj1)$bias, na.rm = TRUE)
    }
    return(obj1)
}

#' @title Reformats raw output from makeChromVAR
#'
#' @description \code{reformatChromVAR} Takes a celltype list of chromVARDeviations or a single chromVARDeviations object over all cell types and reformats it into a ChAI object where each assay is a cell type
#'
#' @param chromVARList The output of makeChromVAR, either as a cell type list (withinCellType and exportList are TRUE, which produces a a named list of ChromVAR objects for each cell type, or a single ChromVAR object run across multiple cell types when withinCellType is FALSE
#' @param selectDev A boolean, with default to FALSE, that decides whether or not deviations or Z score values will be used in the output object. Default of FALSE chooses Z-scores.
#'
#' @return A SummarizedExperiment of either z-scores or deviations across all cell types, (each assay is a cell type), similar to the format of a MOCHA object for scATAC data, or the output of makePseudobulkRNA. 
#'              This format is necessary for using ChromVAR with modeling and association functions (model_General, scATAC_Associations, scRNA_Associations, general_Associations)
#'
#' @keywords data_import
#'
#' @export

reformatChromVAR <- function(chromVARList, selectDev = FALSE){

    CellType <- NULL

    if(!methods::is(chromVARList, 'list') & !methods::is(chromVARList, 'chromVARDeviations')){
        stop('chromVARList is not a list or an individual chromVARDeviation.')
    }else if(methods::is(chromVARList, 'chromVARDeviations')){

        allMeta <-  SummarizedExperiment::colData(chromVARList)
        
        tmpList <- reformatMetaData(allMeta)
        
        sampleData <- tmpList[[1]][sort(rownames(tmpList[[1]])),]
        summarizedData <- tmpList[[2]] 
        
        if(selectDev){
            assayType = 'deviations'
        }else{
            assayType = 'z'
        }
        
        fullMat <- SummarizedExperiment::assays(chromVARList)[[assayType]]
        rownames(summarizedData) <- gsub(" ","_",rownames(summarizedData))
        cellNames <- gsub("__.*","", colnames(fullMat))
       
        assayList <- lapply(unique(cellNames), function(XX){
            
                tmpMat <- fullMat[,cellNames == XX]
                newCellType <- paste(gsub(" |//.","_", XX), collapse ="|")
                colnames(tmpMat) <- sub("__","",gsub(paste("^",newCellType, sep=''),"", colnames(tmpMat)))
                tmpMat[,sort(colnames(tmpMat))]
            
            })
        names(assayList) <- unique(cellNames)
        
    }else if(!all(unlist(lapply(chromVARList, class)) == 'chromVARDeviations')){
        
        stop('Some or all indices of chromVARList are not chromVarDeviations objects. Check your list.')
        
    }else{

        allMeta <- do.call('rbind', lapply(chromVARList, SummarizedExperiment::colData))
        #Identify the end of the original metadata, after which CellType, and various other CellType-Sample metadata was tacked on. 
       
        tmpList <- reformatMetaData(allMeta)

        sampleData <- tmpList[[1]][sort(rownames(tmpList[[1]])),]
        summarizedData <- tmpList[[2]] 

        # Pull out the assays needed and rename them after the celltypes. 
        if(selectDev){
            assayType = 'deviations'
        }else{
            assayType = 'z'
        }
        
        assayList <- lapply(seq_along(chromVARList), function(x){
            tmpMat <- SummarizedExperiment::assays(chromVARList[[x]])[[assayType]]
            newCellType <- paste(gsub(" |//.","_", names(chromVARList)[x]), collapse ="|")
            colnames(tmpMat) <- sub("__","",gsub(paste("^",newCellType, sep=''),"", colnames(tmpMat)))
            tmpMat
            })
        
    }
        
    outputSE <- SummarizedExperiment::SummarizedExperiment(assayList, colData = sampleData, 
                        metadata = list('summarizedData' = summarizedData,
                                        'History' = paste("reformatChromVARList", utils::packageVersion("ChAI"))))

    return(outputSE)
}


#' @title \code{reformatMetaData}
#'
#' @description \code{allMeta} Takes metadata from a combinedSampleTileMatrix, combined ChromVAR object, or combinedPseudobulk object and splits it back out into a sample-specific metadata data.frame and a SummarizedExperiment containing population-specific metrics (celltype by sample). 
#'
#' @param allMeta A data.frame with a CellType and Sample column, containing metadata across all samples and cell type. From the first column untill the Celltype column, sample-specific metadata is encoded, and from the Celltype column till the end, Population-level metadata is encoded.
#'
#' @noRd

reformatMetaData <- function(allMeta){
    
    . <- newCellType <- CellType <- NULL
    
    limitation1 <-which(colnames(allMeta) == 'CellType')
    #Transform metadata for cell type, and find the sample-level data
    newCellType <- paste(gsub(" |//.","_", unique(allMeta$CellType)), collapse ="|")
    newSample <- sub("__","",sub(paste("^",newCellType,sep=''),"", allMeta$Sample))
    originMeta <- allMeta[, 1:(limitation1-1)]
    originMeta$Sample = newSample
    sampleData <- dplyr::distinct(as.data.frame(originMeta))
    rownames(sampleData) <- sampleData$Sample


    ## Transform the CellType-Sample metadata
    sumData <- as.data.frame(allMeta[,(limitation1+1):length(colnames(allMeta))]) 
    sumDatalist <- lapply(colnames(sumData), function(x){
            tmpMat <- sumData[,x, drop = FALSE]
            tmpMat$CellType = allMeta$CellType
            tmpMat$Sample = newSample
            newTmp <- as.data.frame(tidyr::pivot_wider(tmpMat,id_cols = 'CellType', names_from = 'Sample',
                            values_from = {{x}}))
            rownames(newTmp) <- newTmp$CellType
            dplyr::select(newTmp, !CellType)
    })

    names(sumDatalist) <-   colnames(sumData)     

    summarizedData = SummarizedExperiment::SummarizedExperiment(sumDatalist, colData = sampleData)
    return(list(sampleData, summarizedData))

}


#' @title Generate a ChAI scRNA object (Pseudobulk scRNA)
#' 
#' @description \code{makePseuduoulkRNA} pseodubulks a Seurat object by sample and cell type into a SummarizedExperiment, similar to MOCHA

#'
#' @param SO Seurat Object 
#' @param cellTypeColumn The column of the Seurat object with cell type information
#' @param sampleColumn The column of the Seurat object with sample information
#' @param cellPopulations A list of cell type names to extract or the word 'All' (which is default). Determines which cell populations you want to aggregate and analyze from the Seurat object.
#' @param numCores The number of cores to parallelize this operation over. 
#' @param normalize A boolean to determine whether or not to return DESeq2-normalized counts, or raw aggregated counts. Default is TRUE.
#' @param Seurat_format A string, describing the gene name format of the Seurat object. This is used to annotate gene loci, and convert IDs. See documentation for AnnotationDbi::mapIds to identify formats. Default is 'SYMBOL', but 'ENSEMBL' is also common.
#' @param TxDb A Transcript database to be used for identifying gene locations. Must be installed already. Both TxDb and OrgDb must be provided to annote the gene locations.
#' @param OrgDb An Organism database to match up gene names and locations. Must be installed already. Both TxDb and OrgDb must be provided to annote the gene locations.
#' @param verbose Boolean. Default is TRUE, and will print messages. 
#' @return A SummarizedExperiment carrying pseudobulked average expression per 1000 cells for each cell type. 
#'
#' @keywords data_import
#'
#' @export


makePseudobulkRNA <- function(SO, cellTypeColumn, sampleColumn = "Sample", 
                                 cellPopulations = 'All',
                                assayName = 'RNA',
                                numCores = 2,
                                 Seurat_format = 'SYMBOL',
                                 TxDb = NULL, 
                                 OrgDb = NULL,
                                normalize = FALSE,
                             verbose = TRUE) {

    if(any(!c(cellTypeColumn, sampleColumn) %in% colnames(SO@meta.data))){

        stop('Sample or cell type columns are missing from Seurat object. Please verify that the provided column names are correct.')

    }
    
    if(numCores > 3){
    
        warning('User requested to multithread over more than 3 cores. More multithreading will lead to greater memory usage')
    }
    
    if(is.null(TxDb) & !is.null(OrgDb)){

        stop('An OrgDb was provided, but no TxDb was provided. Please provide both, or neither')

    }else if(!is.null(TxDb) & is.null(OrgDb)){

        stop('An TxDb was provided, but no OrgDb was provided. Please provide both, or neither')

    }else if(!methods::is(OrgDb, 'OrgDb') & !is.null(OrgDb)){

        stop('The OrgDb provided is not a Organism Database. Please provide an Organism Database or NULL.')

    }else if(!methods::is(TxDb, 'TxDb') & !is.null(TxDb)){

        stop('The TxDb provided is not a Transcript Database. Please provide a Transcript Database or NULL.')
    }
        
     
    #Generate sample and cell type column in metadata for pseudobulking (after filtering down)
    metadata <- as.data.frame(SO@meta.data)

    ## Detect version of Seurat object
    if(!any(names(SO@assays) %in% assayName)){
        stop('assayName not found within the Seurat object.')
    }else if(methods::is(SO@assays[[assayName]], 'Assay5')){
        counts <- SO@assays[[assayName]]@layers$counts
        cellNames = rownames(SO@assays[[assayName]]@cells)
        colnames(counts) = rownames(SO@assays[[assayName]]@cells)
        rownames(counts) = rownames(SO@assays[[assayName]]@features)
    }else{
        #Seurat v4 or earlier
        counts <- SO@assays[['assayName']]@counts
    }
        
    if(any(dim(counts)==0)){

        stop('Error: Cannot find count matrix within Seurat Object')
    }
    metadata$CellTypeColumn = factor(unlist(metadata[,cellTypeColumn]))
    
    if(any(is.na(SO@meta.data[[sampleColumn]]))){
        stop('Some values within the sampleColumn are NA. Please correct the NAs or choose a different column for sample identification.')
    }
    if(all(rownames(metadata) == c(1:dim(metadata)[1]))){

        stop('Seurat metadata has corrupted rownames. Please make sure the Seurat object has cell names saved.')
        
    }

    if(all(tolower(cellPopulations) == 'all')){
        cellPopulations <- unique(metadata$CellTypeColumn)
    }else if(all(cellPopulations %in% unique(metadata[,cellTypeColumn]))){
        metadata <- dplyr::filter(metadata, !!as.name(cellTypeColumn) %in% cellPopulations)
        counts <- counts[,rownames(metadata)]
    }else{
        stop('Within the cellTypeColumn, not all cellPopulations could be found.')
    }

    cellTypeList <- as.character(unique(metadata$CellTypeColumn))
    sampleList <- unique(SO@meta.data[,sampleColumn])
    emptySample = counts[,1]
    emptySample[TRUE] = 0

    if(verbose){     message("Generating sample List")    }
    
    metadata$Cells = rownames(metadata)
    cellTypeListDFs <- dplyr::group_split(metadata, CellTypeColumn)
    names(cellTypeListDFs) = unlist(lapply(cellTypeListDFs, function(XX) unique(XX$CellTypeColumn)))
    ## Subset down to desired populations
     cellTypeListDFs =  cellTypeListDFs[cellPopulations]
    # Create list to iterate over. 
    subSampleList <- pbapply::pblapply(cl = NULL, X =  cellTypeListDFs, function(XX){
                            subMeta = dplyr::group_split(XX, !!as.name(sampleColumn))
                            subCells = lapply(subMeta, function(ZZ) ZZ$Cells)
                            names(subCells) = unlist(lapply(subMeta, function(ZZ) unique(ZZ[,sampleColumn])))
                            subCells
                        })
    
    if(verbose){  message("Processing total counts and percent detection.")    }
                                                     
    if(numCores > 1){
        cl <- parallel::makeCluster(numCores)
    }else{
        cl = NULL
    }
                                          
    iterList = lapply(subSampleList, function(ZZ){
                    list('Cells' = ZZ, 
                         'Counts' = as.matrix(counts[,which(rownames(metadata) %in% unlist(ZZ))]),
                         'sampleList' = sampleList)
                })
    gc()
    rm(counts)
    countList <- pbapply::pblapply(cl = cl, X = iterList, 
                                   .processCounts, returnCounts = TRUE,
                                   emptySample = emptySample)
    gc()
    percentList <- pbapply::pblapply(cl = cl, X = iterList, .processCounts, 
                                     returnCounts = FALSE, emptySample = emptySample)
                       
    names(countList) <- names(percentList) <- cellTypeList
                         
    if(numCores > 1){
        parallel::stopCluster(cl)
    }
    rm(iterList)
                         
    ##Sort sample-level metadata
    cellColDataNoNA <- BiocGenerics::Filter(function(x) {
        !all(is.na(x))
    }, SO@meta.data)   

    sampleSpecificColumns <- dplyr::group_by(cellColDataNoNA, !!as.name(sampleColumn))  
    sampleSpecificColumns <- dplyr::summarize(sampleSpecificColumns, dplyr::across(dplyr::everything(), dplyr::n_distinct)) 
    sampleSpecificColumns <- dplyr::select(sampleSpecificColumns, tidyselect::where(~all(.x==1)))
    sampleData <- dplyr::distinct(cellColDataNoNA[,c(sampleColumn,colnames(sampleSpecificColumns))])
    # Set sampleIDs as rownames
    rownames(sampleData) <- sampleData[[sampleColumn]]
    sampleData$Sample = sampleData[[sampleColumn]]

    ## Save percent detected for each gene.
    percentDetected <- SummarizedExperiment::SummarizedExperiment(
               percentList,
                colData = sampleData[sampleList,]
    )

    #Get cell numbers
    cellCounts <- as.data.frame(table(metadata[, sampleColumn], metadata[, 'CellTypeColumn']))
    names(cellCounts) <- c(sampleColumn, "CellPop", "CellCount")
    cellCounts <- tidyr::pivot_wider(
        cellCounts,
        id_cols = "CellPop",
        names_from = sampleColumn,
        values_from = "CellCount"
    )
    allCellCounts <- as.data.frame(cellCounts[, -1])
    rownames(allCellCounts) <- cellCounts$CellPop
    allCellCounts <- allCellCounts[sort(cellPopulations),]

    #Process the meta data and extract the total counts and features for each population, as well as cell counts. 
    cellColDataCopy <- data.frame(metadata)

    cellColDataCopy[] <- lapply(cellColDataCopy, function(x) {
        utils::type.convert(as.character(x), as.is = TRUE)
    })

    # Assume all numeric columns are to be saved as additionalCellData
    isNumericCol <- unlist(lapply(cellColDataCopy, function(x) is.numeric(x)))
    additionalCellData <- colnames(cellColDataCopy)[isNumericCol]

    # Group by Sample (rows) and cellPop (columns)
    summarizedData_df <- dplyr::group_by(cellColDataCopy, CellTypeColumn, !!as.name(sampleColumn))
    if (!is.null(additionalCellData)) {
        
        if(verbose){     message("Processing additional metadata.")    }
        
        additionalMetaData <- pbapply::pblapply(cl = NULL, X = additionalCellData, function(x) {
           

            suppressMessages(
                summarizedData2 <- dplyr::summarize(summarizedData_df, meanValues = mean(!!as.name(x), na.rm = TRUE, .groups = "drop")) 
            )

            summarizedData2 <- tidyr::pivot_wider(summarizedData2,
                    id_cols = CellTypeColumn,
                    names_from = !!as.name(sampleColumn),
                    values_from = meanValues
                )

            summarizedData2 <- as.data.frame(summarizedData2)
            rownames(summarizedData2) <- summarizedData2[['CellTypeColumn']]
            summarizedData2 <- summarizedData2[, -1, drop = FALSE]

            # Filter to specific cellPopulations
            summarizedData2 <- summarizedData2[
                rownames(summarizedData2) %in% cellPopulations, ,
                drop = FALSE
            ]

            summarizedData2[sort(cellPopulations),sort(colnames(summarizedData2))]
        })
        names(additionalMetaData) <- additionalCellData

    }else if (is.null(additionalCellData)) {

        additionalMetaData <- NULL

    }
                      
    remove(cellColDataCopy)

    summarizedData <- SummarizedExperiment::SummarizedExperiment(
            append(
            list(
                "CellCounts" = allCellCounts[rownames(additionalMetaData[[1]]),
                                             colnames(additionalMetaData[[1]])]
            ),
            additionalMetaData
            ),
            colData = sampleData[sort(rownames(sampleData)),]
        )

    rnaSE <-  SummarizedExperiment::SummarizedExperiment(countList, colData = sampleData,
         metadata = list('summarizedData' = summarizedData,'detectionRates' = percentDetected,
            History = paste("makePseudobulkRNA", utils::packageVersion("ChAI"))))

    #Decide if you are interweaving this data with genomic databases/locations. 
                                  
    #Pull in Transcript and Organism databases. 
    if(!is.null(TxDb) & !is.null(OrgDb)){
        
        newSE <- linkToGenome(rnaSE = rnaSE, gene_format = Seurat_format, 
                              TxDb= TxDb, OrgDb = OrgDb)
       
    }else{
        newSE <- rnaSE

    }

    if(normalize){
        newSE <- normalizePseudobulk(newSE, sampleColumn = sampleColumn)
    }
    return(newSE)
}
                                  
                                  
#' @title \code{subSampleList}
#'
#' @description \code{subSampleList} Helper function for makePseudobulkRNA. 

#'
#' @param cellType  the output of makePseudobulkRNA a
#' @param sampleList1 List of sample names
#' @param sampleColumn1 Name of metadata column with sample names
#' @param metadata1 Seurat metadata
#' @return a list of metadata by celltype and sample
#'
#' @noRd

.subSampleList <- function(cellType, sampleList1, sampleColumn1, metadata1){
        
        cellNameList <- lapply(sampleList1, function(z){
            rownames(dplyr::filter(metadata1, 
                            CellTypeColumn == cellType & !!as.name(sampleColumn1) == z))
        })
        return(cellNameList)
}

#' @title \code{processCounts}
#'
#' @description \code{processCounts} Helper function for makePseudobulkRNA. 

#'
#' @param sampleList_index One index from the subSampleList
#' @param returnCounts Boolean. If true, will return count matrix. If false, will return percent detected. 
#' @param emptySamples A matrix representing an empty sample (all genes 0)
#' @return a list of matrices, either counts or percent detected
#'
#' @noRd
                   

.processCounts <- function(sampleList_index, returnCounts = TRUE, emptySamples){
    
    if(returnCounts){     

    matList <- lapply(sampleList_index[[1]], function(YY){
                        
            if(length(YY) > 1){
                    unlist(rowSums(sampleList_index[[2]][,YY]))
            }else if(length(YY) == 1){
                sampleList_index[[2]][,YY]
            }else{
                emptySample
            }

        }) 

    }else{

        matList <- lapply(sampleList_index[[1]], function(YY){

                if(length(YY) > 1){
                        unlist(rowSums(sampleList_index[[2]][,YY] > 0))/length(YY)
                }else if(length(YY) == 1){
                    as.integer(sampleList_index[[2]][,YY] > 0)
                }else{
                    emptySample
                }

            }) 
	
    }

    mat <- do.call('cbind', matList)
    colnames(mat) <- names(sampleList_index[[1]])
    ## Identify and fill in missing samples with NAs. 
    if(any(!sampleList_index[[3]] %in% colnames(mat))){
        
        missingSamples = sampleList_index[[3]][!sampleList_index[[3]] %in% colnames(mat)]
        
        for(newCol in missingSamples){
        
            mat = cbind(mat, newCol = 0)
            colnames(mat)[colnames(mat) == 'newCol'] = newCol
        }
        
    }
    mat = mat[,sampleList_index[[3]]]
    
    return(mat)
}

#' @title \code{linkToGenome}
#'
#' @description \code{linkToGenome} Helper function for makePseudobulkRNA. Links genes to genomic loci, and tosses poorly annotated transcripts. 
#'
#' @param rnaSE the output of makePseudobulkRNA
#' @param gene_format A string, describing the gene name format of the Seurat object. This is used to annotate gene loci, and convert IDs. See documentation for AnnotationDbi::mapIds to identify formats. Default is 'SYMBOL', but 'ENSEMBL' is also common.
#' @param TxDb A Transcript database to be used for identifying gene locations. Must be installed already. Both TxDb and OrgDb must be provided to annote the gene locations.
#' @param OrgDb An Organism database to match up gene names and locations. Must be installed already. Both TxDb and OrgDb must be provided to annote the gene locations.
#' @return a list of metadata by celltype and sample
#'
#' @keywords data_import
#'
#' @export

linkToGenome <- function(rnaSE, TxDb= NULL, OrgDb = NULL, gene_format = 'SYMBOL'){

    ## fix global bindings
    ALIAS <- SYMBOL <- NULL
    if(!is.null(TxDb) & !is.null(OrgDb)){
        if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
                stop(
                "Package 'AnnotationDbi' is required for adding gene location information.",
                "Please install AnnotationDbi'from Bioconductor to proceed."
                )
        }

        #Extact the GenomicRanges for all genes, while filtering out non-standard chromosomes that mess things up. 
        txList <- suppressMessages(suppressWarnings(S4Vectors::stack(GenomicFeatures::genes(TxDb,
                                                single.strand.genes.only = FALSE))))
        txList <- suppressMessages(GenomeInfoDb::keepStandardChromosomes(sort(txList), 
                                    species='Homo_sapiens',
                                        pruning.mode = 'coarse'))
        #Reduce ranges merges transcripts for the same gene into one, so that we can one ball-park stop and end. 
        txList <- plyranges::reduce_ranges_directed(plyranges::group_by(txList, gene_id))
        txList$GeneSymbol <- suppressWarnings(AnnotationDbi::mapIds(OrgDb, as.character(txList$gene_id), gene_format, "ENTREZID"))
        txList <- plyranges::filter(txList, !is.na(GeneSymbol))
        txList$MultipleMappings = txList$GeneSymbol %in% txList$GeneSymbol[duplicated(txList$GeneSymbol)]
        txList <- plyranges::slice(plyranges::group_by(txList, GeneSymbol), 1)
        names(txList) <- txList$GeneSymbol
        txList <- plyranges::ungroup(txList)
        promoters <- GenomicRanges::promoters(txList)
        promoters <- paste(GenomicRanges::seqnames(promoters), ":", GenomicRanges::start(promoters),"-", GenomicRanges::end(promoters), sep = '')
        txList$Promoters = promoters
        txList$Strand_Direction = as.character(GenomicRanges::strand(txList))
        ## Some genes are mapped to both positive and negative strands, and some miRNAs are provisionally mapped to multiple locations. 
        ## Let's just take the first location for these. 
        
        if(any(!rownames(rnaSE)  %in% txList$GeneSymbol)){

            ## if it's a gene symbol, sometimes there's mismatches, with multiple potential gene names. 
            ## let's make sure that isn't the case, and match up names that don't match the database for SYMBOL. 
            if(gene_format == 'SYMBOL'){
               lostTxs <- rownames(rnaSE)[!rownames(rnaSE)  %in% txList$GeneSymbol]
                #Check for missed due to '.2' gene name versions
               
               aliasAnnot <- AnnotationDbi::select(OrgDb, keys=lostTxs,
                        columns=c("SYMBOL","ENTREZID"), 
                                                   keytype="ALIAS")
                # Filter to matching SYMBOL
                # Identify all ALIASes that were matched well
                aliasAnnot1  <- dplyr::filter(aliasAnnot , !is.na(SYMBOL))
                # Identify all ALIASes that didn't. 
                # Remove .[0-9] and try again.
                aliasAnnot2  <- dplyr::filter(aliasAnnot , is.na(SYMBOL))
                aliasAnnot2  <- dplyr::mutate(aliasAnnot2, ALIAS2 =
                                              gsub("\\.[0-9].*","", ALIAS))
                aliasAnnot3 <- AnnotationDbi::select(OrgDb, 
                            keys=unique(aliasAnnot2$ALIAS2),
                        columns=c("SYMBOL","ENTREZID"), 
                                                   keytype="ALIAS") 
                aliasAnnot3 <- dplyr::filter(aliasAnnot3, !is.na(SYMBOL))
                aliasAnnot3 <- dplyr::inner_join(
                                    aliasAnnot2[,c('ALIAS', 'ALIAS2')],
                                    aliasAnnot3, 
                            by = c('ALIAS2' = 'ALIAS'))[,
                                        c('ALIAS','SYMBOL','ENTREZID')]
                aliasAnnot4 <- dplyr::filter(rbind(aliasAnnot1, aliasAnnot3),
                                             SYMBOL  %in% txList$GeneSymbol &
                                            !SYMBOL %in% rownames(rnaSE))
                aliasAnnot5 <- dplyr::slice_head(
                                    dplyr::group_by(aliasAnnot4, ALIAS),
                                    n=1)
              
                subTx <- txList[txList$GeneSymbol %in% aliasAnnot5$SYMBOL]
                subTx$Alias = subTx$GeneSymbol
                subTx$GeneSymbol = aliasAnnot4$ALIAS[match(subTx$GeneSymbol, 
                                                       aliasAnnot5$SYMBOL)]
                names(subTx) = subTx$GeneSymbol
                subTx$Alias = aliasAnnot4$ALIAS[match(subTx$GeneSymbol, 
                                                       aliasAnnot5$SYMBOL)]
                txList = sort(c(txList, subTx))

              
            }
            #Identify Transcripts that were removed.
            lostTxs <- rownames(rnaSE)[!rownames(rnaSE)  %in% txList$GeneSymbol]
            
            #Subset down the read count matrix to just the transcripts we have in this database.
            newSE <- rnaSE[rownames(rnaSE)  %in% txList$GeneSymbol, ]
            
            message(stringr::str_interp("${length(lostTxs)} transcripts were not found within the provided transcript/organism database, or               did not have a clear genomic location within the standard chromosome assembly for your organism. 
                Filtered transcript IDs are saved in metadata, via obj@metadata[['filteredIDs']]."))
                    
            newSE@metadata = append(newSE@metadata, list('filteredIDs' = lostTxs))

        }
            
        #Subset down the read count matrix to just the transcripts we have in this database.
        newSE <- rnaSE[rownames(rnaSE)  %in% txList$GeneSymbol, ]
        
        # Repackage the gene-sample matrices, the sampleData, transcript Genomic Ranges, and associated metadata (countInfo) into one SummarizedExperiment object. 
        attachedList <- plyranges::ungroup(txList)
        SummarizedExperiment::rowRanges(newSE) = attachedList[match(rownames(newSE), names(attachedList))]
       
    }else{
        newSE <- rnaSE
        warning('No TxDb and/or OrgDb provided.')
    }
    return(newSE)
}



#' @title Normalize a ChAI scRNA object using DESEq2
#'
#' @description \code{normalizePseudobulk} Takes the output of makePseudobulkRNA and normalizes it. 

#'
#' @param rnaSE  the output of makePseudobulkRNA a
#' @param sampleColumn The column of the Seurat object with sample information
#' @return A SummarizedExperiment with normalized average expression
#'
#' @keywords data_import
#'
#' @export

normalizePseudobulk <- function(rnaSE, sampleColumn = 'Sample'){

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop(
        "Package 'DESeq2' is required for normalizePseudobulk. ",
        "Please install 'DESeq2' to proceed."
        )
    }

    allMat <- lapply(names(SummarizedExperiment::assays(rnaSE)), function(x){
                print(x)
                old_mat <- SummarizedExperiment::assays(rnaSE)[[x]]
                old_mat[is.na(old_mat)] = 0
                old_mat = round(old_mat) #ensure integers
                suppressMessages(
                    dds <- DESeq2::DESeqDataSetFromMatrix(countData = old_mat[,colSums(old_mat) > 0],
                              colData = SummarizedExperiment::colData(rnaSE)[colSums(old_mat) > 0,],
                                design = stats::as.formula(paste( '~', sampleColumn)))
                )

                suppressMessages(dds <- DESeq2::estimateSizeFactors(dds))
                suppressMessages(new_mat <- DESeq2::counts(dds, normalize = TRUE))
                if(any(colSums(old_mat) == 0)){
                    filled_data <- do.call('cbind', lapply(which(colSums(old_mat) == 0), function(x){
                            rep(0, dim(new_mat)[1])
                    }))
                    new_mat <- cbind(new_mat, filled_data)
                    
                }
                
                new_mat <- new_mat[,sort(colnames(new_mat))]
        })

    names(allMat) <- names(SummarizedExperiment::assays(rnaSE))
    
    se <- SummarizedExperiment::SummarizedExperiment(allMat, 
                        colData =SummarizedExperiment::colData(rnaSE)[sort(colnames(rnaSE)),],
                        metadata = rnaSE@metadata)
    SummarizedExperiment::rowRanges(se) <-  SummarizedExperiment::rowRanges(rnaSE)

    se@metadata$History <- append(se@metadata$History, paste("normalizePseudobulk", utils::packageVersion("ChAI")))
    return(se)
    
    }

CellTypeColumn <- gene_id <- GeneSymbol <- meanValues <- NULL





#' @title Functions for identifying genes that pass threshold for a given cell types. 
#'
#' @description \code{thresholdGenes} Identifies which genes pass the threshold settings for a given cell type
#' @param rnaSE A scRNA SummarizedExperiment from makePseudobulkRNA. 
#' @param factors Categorical factors, found within the sample metadata. If numeric factors are provided,
#'                   they will be avoided. If genes pass the provided threshold in any group, 
#'                   they will be returned.
#' @param detectionThreshold A number between 0 and 1, representing the mean detection rate threshold for a given gene to be modeled. 
#'   This detection rate is calculated for each sample and cell type during makePseudobulkRNA, and represents the percentage of cells that have a transcript for a given gene. Over all samples, the average has to be above this to be modeled. Default is 0.01. 
#' @param expressionThreshold A number greater than zero, representing the expression threshold for modeling. A given gene, on average across all samples, expressed above this threshold. The default is 0. 
#' @param cellCountThreshold The minimum number of cells in a given pseudobulk for the pseudobulk to be included in analysis. If fewer than this number of cells are found, then the sample will be dicarded The number of cells within the pseudobulked scRNA. Default is 10 cells. 
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return output_vector A linear model
#'
#' @export
#' @keywords data_object_manipulation

thresholdGenes <- function(rnaSE,   
                           factors = NULL,
                        cellPopulation = 'ALL',
                    detectionThreshold = 0.01,
                    expressionThreshold = 0,
                    cellCountThreshold = 10){
    
  if(isChAIObject(rnaSE, type = 'data', returnType = TRUE) != 'scRNA'){

    stop("rnaSE is not a ChAI scRNA Object (i.e. normalized, pseudobulked via ChAI's makePseudobulkRNA function)")

  }

    # check cell Populations
  if(!all(cellPopulation %in% names(SummarizedExperiment::assays(rnaSE))) &
     all(tolower(cellPopulation) != 'all')){
      
    stop('rnaSE does not contain an assay that matches the cellPopulation input variable.')
      
  }else if(all(tolower(cellPopulation) == 'all')){
    
          cellPopulation = SummarizedExperiment::assayNames(rnaSE)
  }

 if(all(cellPopulation != 'counts')){
      detectionRates = rnaSE@metadata$detectionRates
      summarizedData = rnaSE@metadata$summarizedData
  }else{
      detectionRates = rnaSE@metadata$detectionRates
      summarizedData = SummarizedExperiment::colData(rnaSE)
    }

      
  if(!methods::is(cellCountThreshold, 'numeric')){
 
      stop('cellCountThreshold must be numeric.')
      
 }

  metadf = SummarizedExperiment::colData(rnaSE)
      
  if(!is.null(factors) & any(factors %in% colnames(metadf))){

      charClass = factors[unlist(lapply(factors, function(XX) { 
          class(metadf[,XX]) %in% c('character','factor')}))]
      
 }else{
    
      charClass = character(0)
      
    }
    
  geneList = lapply(cellPopulation, function(cellPop){
          if(cellPop != 'counts'){
              passThreshold = SummarizedExperiment::assays(summarizedData)[['CellCounts']][cellPop,] > cellCountThreshold 
          }else{
              passThreshold = summarizedData$CellCounts > cellCountThreshold 
          }
          subDetect = detectionRates[passThreshold,]
          subRNA = rnaSE[passThreshold,]
      
          if(length(charClass) == 0){
      
              detectMean <- rowMeans(SummarizedExperiment::assays(subDetect)[[cellPop]])
              expressMean <- rowMeans(SummarizedExperiment::assays(subRNA)[[cellPop]])
              
              intersect(names(which(expressMean > expressionThreshold)),
                                              names(which(detectMean > detectionThreshold)))
              
            }else{

                 #Identify groups for a given categorical variable. 
              unique(unlist(lapply(charClass, function(XX){
                    groups = unique(metadf[passThreshold,XX])
                    groups = groups[!is.na(groups)]
                    unique(unlist(lapply(groups, function(ZZ){
                                   #Test to see which genes pass threshold within that group
                                    specificSamples = metadf[,XX] == ZZ
                                    specificSamples[is.na(specificSamples)] = FALSE
                                    detectMean <- rowMeans(SummarizedExperiment::assays(
                                        subDetect[,specificSamples])[[cellPop]])
                        
                                    expressMean <- rowMeans(SummarizedExperiment::assays(
                                        subRNA[,specificSamples])[[cellPop]])
                        
                                    intersect(names(which(expressMean > expressionThreshold)),
                                              names(which(detectMean > detectionThreshold)))
                              })))
                })))
            }
            
      })
     names(geneList) = cellPopulation
    return(geneList)
}
                                         
isChAIObject <- function(Object, type = 'data', returnType = FALSE) {
   
  if (methods::is(Object, "SummarizedExperiment")){

    if(!any(names(Object@metadata) %in% 'History') & type == 'data'){
         stop("Object was not generated by ChAI or MOCHA.")
    } else if(!any(grepl('getSampleTileMatrix|makePseudobulkRNA|importGeneralModality|reformatChromVARList', unlist(Object@metadata$History))) & type == 'data'){
      stop("Object is not an SampleTile object from MOCHA, or a SummarizedExperiment object from ChAI's makePseudobulkRNA or importGeneralModality.")
    }else if(!any(grepl('getSampleTileMatrix|makePseudobulkRNA|importGeneralModality|reformatChromVARList', unlist(Object@metadata$History))) & type == 'data'){
      stop("Object is not an SampleTile object from MOCHA, or a SummarizedExperiment object from ChAI's makePseudobulkRNA, makeChromVAR, or importGeneralModality.")
    }else if(! type  %in% c('data', 'model')){
      stop("type not recognized. Must be either 'data' or 'model'")
    }else if(type == 'data' & returnType){

      specObject <- unlist(lapply(c("transformChAI", 'getSampleTileMatrix',
                    'makePseudobulkRNA','importGeneralModality', 'reformatChromVARList', 
                    'runSSGSEA'), function(XX){
        
          any(grepl(XX, Object@metadata$History)) 
        }))

      if(specObject[6]){
        return('General')
      }else if(specObject[1]){
        return('Transformed')
      }else if(specObject[2]){
        return('scATAC')
      }else if(specObject[3]){
        return('scRNA')
      }else if(specObject[4]){
        return('General')
      }else if(specObject[5]){
        return('ChromVAR')
      }else{
        return('model')
      }
            
    }else if(type == 'data'){

          specObject <- unlist(lapply(c("transformChAI", 'getSampleTileMatrix',
                    'makePseudobulkRNA','importGeneralModality', 'reformatChromVARList'), function(XX){
        
          any(grepl(XX, Object@metadata$History)) 
        }))

        if(any(specObject)){

            return(TRUE)
            
        }

    }

    if(!any(names(Object@metadata) %in% 'Type') & type == 'model'){
         stop("Object is not a model object generated by ChAI.")
    } else if(!any(grepl(paste(getModelTypes(), collapse = '|'), Object@metadata$Type)) & type == 'model'){
         stop("Object is not a ChAI model object.")
    }else if(type == 'model' & returnType){
      return(Object@metadata$Type)
    }
    return(FALSE)

  }else if(methods::is(Object, "list")){
        
       allModels = lapply(Object, function(XX) isChAIObject(XX, type = type, returnType = returnType))
       allModels = unique(unlist(allModels))
       allModels = allModels[allModels != 'FALSE']
       if(length(allModels) > 1){ return(FALSE)}
       return(allModels)
          
    }else{
      return(FALSE)
  }
}                                                                           