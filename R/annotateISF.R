#' Function to annotate features using LC-MS in-source fragmentation 
#' (ISF) information.
#'
#' Search possible matches of a feature in Lipid fragments and other small
#' molecule libraries, using 'raw' LC-MS single chromatogram information
#' already loaded and peak-picked.
#'
#' @author Goncalo Graca (Imperial College London)
#'
#' @param targetTable A .csv file containing the list of features to annotate
#' and the name of the files containing the raw data.
#' @param filetype LC-MS chromatogram file extension ("mzML" or "CDF"). 
#' The default is "mzML".
#' @param ESImode Ionization mode: 'POS' for positive (default)
#' or 'NEG' for negative ionisation modes
#' @param libs Fragment libraries to use: 'Lipids' (default) or 'Metabolites'
#' for other small molecules.
#' @param RTfile Optional csv file with Lipid/metabolites classes Retention
#' Times in seconds.
#' @param corThresh Pearson correlation coefficient for EIC correlation.
#' @param checkIsotope Whether or not to check the isotope type;
#' default is set to TRUE
#' @param tolerance Tolerance in ppm for the candidate search.
#' @param maxMZdiff Maximum m/z difference between candidate fragments and
#' pseudo-MS/MS or MS scans ions in Da.
#' @param matchWeight weight of the fragment matches to the final score;
#' value between 0 and 1; the remaining fraction of the weight comes from the
#' candidate m/z error.
#' @param ncandidates Maximum number of candidates to plot and store.
#' @return For each feature in the targetTable the will return a ranked list of
#' annotations, a plot of EICs and ISF pseudo-MS/MS spectrum for the matched 
#' ions, a plot of pseudo-MS/MS and pseudo-MS spectra for the each feature,
#' a targetTable annotated with rank 1 annotations and a table with the options
#' used for the function.
#' @examples
#' # Get the some example human serum LC-MS data and feature list to annotate:
#' getDemoData()
#' # Download the example .mzML from zenodo website and store it in the 
#' # working directory:
#' # https://zenodo.org/records/17408169
#' # Run the annotation using the lipid libraries:
#' annotateISF(targetTable="targetTable.csv", filetype="mzML", 
#' libs="Lipids", ESImode="POS", RTfile="none", corThresh=0.8,
#' checkIsotope=TRUE)
#' @importFrom utils read.csv write.csv
#' @importFrom grDevices dev.off pdf
#' @export
annotateISF <- function(targetTable="targetTable.csv",
                        filetype="mzML",
                        libs="Lipids",
                        ESImode="POS",
                        RTfile="none",
                        corThresh=0.8,
                        checkIsotope=TRUE,
                        tolerance=25,
                        maxMZdiff=0.01,
                        matchWeight=0.5,
                        ncandidates=5){
    ## Read XCMS peak-picking options------------------------------------------
    xcmsOptions <- getXcmsOptions()
    if(is.null(xcmsOptions)){
        stop("Please edit the 'xcmsOptions.csv' and re-run.")
    }
    ## Read targets table------------------------------------------------------
    if(!file.exists("targetTable.csv")) {
        stop("targetTable.csv not found! Run getDemoData(), 
                edit targetTable.csv and re-run.")
    }
    targets <- read.csv("targetTable.csv")
    
    ## Initialize results object-----------------------------------------------
    results <- initializeResultsISF(targets, libs, ESImode)
    
    ## Get RTs intervals from file --------------------------------------------
    RTs <- RTsFromFile(RTfile)
    
    ## load libraries ---------------------------------------------------------
    libraries <- loadLibs(libs, ESImode)
    libfiles <- libraries$libfiles
    lib <- libraries$lib
    
    # evaluate if one or more samples are to be read---------------------------
    if(length(unique(targets[,3])) == 1){
        message("Reading data...")
        MSData <- readMS1(filetype=filetype, target=targets[1,3], 
                            xcmsOptions=xcmsOptions)
        xcmsF1 <- MSData$xcmsF1
        peaksF1 <- MSData$peaksF1
    } else NULL

    ## process each feature from the targets table-----------------------------
    for(i in seq_len(nrow(targets))){
        progNote <- paste("####### Processing feature", i, "of", nrow(targets),
                            "#######")
        message(progNote)
        if(length(unique(targets[,3]))>1){
            message("Reading data...")
            MSData <- readMS1(filetype=filetype, target=targets[i,3], 
                            xcmsOptions=xcmsOptions)
            xcmsF1 <- MSData$xcmsF1
            peaksF1 <- MSData$peaksF1
        } else NULL
    
        fmz <- targets[i,1]
        frt <- targets[i,2]
    
        ## get sample name to save files into --------------------------------
        SpName <- targets[i,3]
    
        ## get MS spectra at feature RT --------------------------------------
        message("Obtaining ISF spectrum...")
        try(
            specs <- getPseudoMSMS(fmz, frt, xcmsF1, xcmsF2=NULL, 
                            peaksF1, peaksF2=NULL,filetype=filetype, 
                            nCE=1, cthres1=corThresh,
                            cthres2=corThresh, savePlotResults=TRUE, 
                            savePseudoMSMS=TRUE, ExpName="LCMS", 
                            DirPath=results$resultsDir)
        )
    
        if(exists("specs")){
            highCESpec <- specs$ms1
            pseudoSpec <- specs$insource
            inSourceSpec <- specs$insource
            isfEIC <- specs$ms1_eic
            feic <- specs$feic
        } else next
    
        ## Isotope check ------------------------------------------------------
        if(!checkIsotope) iso <- 0 else {
            #message("Checking isotope type...")
            if(is.null(inSourceSpec)) {
                iso <- NA
        } else if(length(inSourceSpec) <= 12) {
                iso <- 0
        } else {
                iso <- checkIsotope(fmz, frt, inSourceSpec)
            }
        }
    
        ## Search Libraries ---------------------------------------------------
        if(is.null(pseudoSpec) & is.null(highCESpec)) { next
            } else {
            message("Searching candidates...")
            candidates <- searchLib(lib, libfiles, fmz - iso * 1.0034,
                                frt, tolerance = tolerance, RTs, inSourceSpec)
        }
    
        # Compare fragments between Library candidates and
        # high-collision-energy / pseudo-MS/MS spectra ------------------------
        if(is.null(pseudoSpec) & is.null(highCESpec) |
            length(unlist(candidates)) == 0) {
            result <- NULL
        } else {
        message("Matching candidate(s) fragments to ISF spectra...")
        output <- mapply(
            compFrag, candidates, lapply(as.numeric(names(candidates)), 
            function(x) lib[[x]]), 
            MoreArgs=list(fmz, frt, iso, highCESpec, pseudoSpec, 
                            maxMZdiff=maxMZdiff, matchWeight=matchWeight), 
            SIMPLIFY = FALSE)
        result <- do.call(rbind,lapply(output, "[[", 1))
        specMatch <- unlist(lapply(output, "[[",2), recursive = FALSE)
        specMatch <- specMatch[!(specMatch) == "NULL"]
        }
        ## Score ranking ------------------------------------------------------
        if(is.null(result)) {
        rankedResult <- targets[i, c(1,2)]
        rankedResult[c("metabolite", "feature.type", "ion.type", "isotope",
                    "mz.metabolite", "matched.mz", "mz.error", "pseudoMSMS",
                    "fraction", "score")] <- NA
        # type of ion isotope
        rankedResult$isotope <- paste("M+", iso, sep="")
    
        # pseudoMSMS flag
        if(is.null(pseudoSpec)) {
                rankedResult$pseudoMSMS <- "FALSE"
            } else { rankedResult$pseudoMSMS <- "TRUE"
            }
        } else {
            output <- rankScore(result,specMatch)
            rankedResult <- output$rankedResult
            rankedSpec <- output$rankedSpecMatch
        }

        ## Save output of ranked annotations ----------------------------------
        if(exists("rankedResult")){
            saveRankedISF(fmz, frt, ncandidates, rankedResult, SpName, 
                    resultsDir=results$resultsDir)
        }

        ## Store highest rank annotation in global results---------------------
        if(!exists("rankedResult")) rankedResult <- NULL
        results$global[i,] <- storeAnnotations(global=results$global[i,], 
                                                rankedResult[1,])

        ## save plot of matched spectra for the top n candidates---------------
        if(exists("rankedSpec")){
            saveMatchedISF(fmz, frt, highCESpec, result, isfEIC, output, 
                            ncandidates, rankedSpec, SpName, 
                            resultsDir=results$resultsDir)
        }
    }
    ## save global results table-----------------------------------------------
    write.csv(results$global,
                file = paste(results$resultsDir, "Global_Results",
                ".csv", sep = ""), row.names = FALSE)
    # save general options
    df <- data.frame(targetsTable_file=targetTable, libraries=libs,
                        ESImode=ESImode, RTfile=RTfile, corThresh=corThresh,
                        checkIsotope=checkIsotope, matchWeight=matchWeight,
                        tolerance=paste(tolerance, "ppm"),
                        maxMZdiff=paste(maxMZdiff, "Da"), row.names="Option")
    df <- as.data.frame(t(df))
    write.csv(df, file=paste(results$resultsDir, "General_options", 
                                ".csv", sep = ""))
    # save xcms options
    write.csv(xcmsOptions,
                file=paste(results$resultsDir, "XCMS_options", ".csv", sep=""), 
                row.names=FALSE)
    message("Job done!")
}

## Helper functions------------------------------------------------------------

## Initialization of results folder and global results table-------------------
initializeResultsISF <- function(targets, libs, ESImode){
    # Create directory to store the results
    mainDir <- "./Annotations"
    Date <- Sys.Date()
    Time <- format(Sys.time(), "%X")
    Time <- gsub(":", "_", Time)
    subDir <- paste(libs, "_", ESImode, "_ISF_", Date,"_", Time, sep="")
    dir.create(file.path(mainDir), showWarnings = FALSE)
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    resultsDir <- paste(mainDir, "/", subDir, "/", sep = "")
    # create table to store global results
    global <- targets
    global[,c("metabolite", "feature.type", "ion.type", "isotope", 
                "mz.metabolite","matched.mz", "mz.error", "pseudoMSMS", 
                "fraction", "score")] <- NA
    # return global results table and results path as list
    results <- list(global=global, resultsDir=resultsDir)
    return(results)
}

## Load table with XCMS peak-picking options-----------------------------------
getXcmsOptions <- function(){
    if(file.exists("XCMS_options.csv")) {
        xcmsOptions <- read.csv("XCMS_options.csv")
    } else {
        message("No XCMS options file found in the working directory")
        answer <- readline(prompt = "Use default options (Yes/No)? [y/n]")
        if(answer == "y") {
            xcmsOptionsPath <- system.file("XCMS_options.csv",
                                            package = "MetaboAnnotatoR")
            xcmsOptions <- read.csv(xcmsOptionsPath)
        }
        if(answer == "n"){
            xcmsOptionsPath <- system.file("XCMS_options.csv",
                                            package = "MetaboAnnotatoR")
            file.copy(from = xcmsOptionsPath, to = getwd())
            xcmsOptions <- NULL
            message("Default  XCMS options file saved in the working directory"
                    )
        }
    }
    return(xcmsOptions)
}

## read data from mzML or CDF files--------------------------------------------
readMS1 <- function(filetype, target, xcmsOptions){
    if(filetype == "mzML"){
        dataPath <- paste(target, ".mzML", sep = "")
        # separate the two MS functions
        xcmsF1 <- MSnbase::readMSData(dataPath, msLevel. = 1, mode = "onDisk")
    }
    if(filetype == "CDF"){
        dataPath <- paste("./", targets[1,3], "01.CDF", sep = "")
        xcmsF1 <- MSnbase::readMSData(dataPath, mode = "onDisk")
    }
    # centwave parameters
    cwp <- xcms::CentWaveParam(snthresh = xcmsOptions[3,2],
                                noise = xcmsOptions[2,2],
                                ppm = xcmsOptions[1,2],
                                peakwidth = as.numeric(xcmsOptions[4,2:3]),
                                prefilter = as.numeric(xcmsOptions[5,2:3]))
    peaksF1 <- xcms::findChromPeaks(xcmsF1, msLevel = 1L, param = cwp)
    # gather all data objects into a list
    MSData <- list(xcmsF1=xcmsF1, peaksF1=peaksF1)
    return(MSData)
}

## Save output of ranked annotations ------------------------------------------
saveRankedISF <- function(fmz, frt, ncandidates, rankedResult, 
                            SpName, resultsDir){
    write.csv(rankedResult[rankedResult$rank <= ncandidates,],
                file=paste(resultsDir, SpName, "_", round(fmz,3), "mz_", 
                        round(frt,3), "s_", "ranked_candidates.csv", 
                        sep=""), 
                row.names=FALSE)
}

## save EICs and pseudoMSMS spectra of matched candidates----------------------
# save image with matched spectra and EICs
# must update the folder to save images in...
# can be updated to work with single feature
saveMatchedISF <- function(fmz, frt, 
                        highCESpec,
                        result,
                        isfEIC,
                        output, 
                        ncandidates,
                        rankedSpec,
                        SpName,
                        resultsDir){
    if(!is.null(result)){
        if(length(rankedSpec) == 1) {
            plotCandidatesAIF(fmz, frt, highCESpec, 
                                isfEIC, SpName, output, 1, resultsDir)
        }
        if(length(rankedSpec) <= ncandidates) {
            plots <- lapply(seq_len(length(rankedSpec)),
                            function(x) plotCandidatesAIF(fmz, frt,
                                                            highCESpec,
                                                            ms2eic=isfEIC, 
                                                            SpName,
                                                            output, x, 
                                                            resultsDir))
        }
        if(length(rankedSpec) > ncandidates) {
            plots <- lapply(seq_len(length(rankedSpec[seq_len(ncandidates)])),
                            function(x) plotCandidatesAIF(fmz, frt,
                                                            highCESpec,
                                                            ms2eic=isfEIC, 
                                                            SpName,
                                                            output, x, 
                                                            resultsDir))
        }
    } else NULL
}
