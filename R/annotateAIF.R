#' @title Annotate features from LC-MS AIF raw chromatograms
#'
#' @description
#' This function annotates features from raw LC-MS AIF chromatograms, by 
#' performing pseudo-MS/MS spectra deconvolution and then matching ions 
#' to metabolite/lipid fragment libraries.
#'
#' @author Goncalo Graca (Imperial College London)
#'
#' @param targets A data frame containing the features to annotate
#' and the file paths to the raw data.
#' @param xcmsOptions A data frame containing the XCMS \code{centWave} 
#' peak-picking parameters. An example of such table can be found in the data
#' provided with \code{MetaboAnnotatoR} as \code{XCMS_options.csv} 
#' (see example for details).
#' @param libs Fragment libraries to use. Either the built-in libraries can be
#' specified (\code{LipidPOS}, \code{LipidNEG}, \code{MetabolitesPOS}, 
#' \code{MetabolitesNEG}) or the full path to user-defined libraries.
#' @param RTs Optional data.frame with Lipid/metabolites classes Retention
#' Times in seconds.
#' @param nCE Number of Collision Energy levels depending on the MS system used
#' Waters, Bruker (QToF) and Thermo Orbitrap = 1, Agilent (QToF) > 1, however,
#' only the highest energy level will be considered.
#' @param corThresh Pearson correlation coefficient for EIC correlation.
#' @param checkIsotope Whether or not to check the isotope type;
#' default is set to TRUE
#' @param tolerance Tolerance in ppm for the candidate search.
#' @param maxMZdiff Maximum m/z difference between candidate fragments and
#' pseudo-MS/MS or AIF ions in Da.
#' @param matchWeight weight of the fragment matches to the final score;
#' value between 0 and 1; the remaining fraction of the weight comes from the
#' candidate m/z error.
#' @return For each feature in the targets data table the function will return 
#' a data frame with each feature rank 1 annotation and a table with the 
#' options used for the function, the data and time of annotation. In addition, 
#' lists with the ranked candidates matched to each feature
#' (\code{rankedResult}), ranked matched spectra (\code{rankedSpectra}) and a 
#' list with pseudo-MS/MS spectra, in-source spectra, and AIF spectra and 
#' respective EIC objects (\code{pseudoMSMS}, see also \code{getPseudoMSMS},
#' documentation for details).
#' @examples
#' # Get the example feature list and peak-picking parameters files:
#' # Download the example .mzML from zenodo website into the working directory:
#' download.file(
#' "https://zenodo.org/records/17408169/files/Lipid_Positive_QC.mzML?download=1", 
#' "Lipid_Positive_QC.mzML"
#' )
#' # create a new targetTable with one feature to annotate
#' targets <- data.frame(feature.mz=520.3408533, feature.rt=100.6238759, 
#' Sample.name="Lipid_Positive_QC.mzML")
#' # read the default xcms parameters on the XCMS_options.csv file and modify
#' # the noise threshold parameter
#' xcmsOptionsPath <- system.file("XCMS_options.csv", 
#' package="MetaboAnnotatoR")
#' xcmsOptions <- read.csv(xcmsOptionsPath)
#' xcmsOptions[2,2] <- 1000
#' # Run the annotation using the built-in lipid POS library:
#' annotations <- annotateAIF(targets, xcmsOptions, 
#' libs="LipidPos", RTs="none", nCE=1, corThresh=0.8,
#' checkIsotope=TRUE)
#' @export
annotateAIF <- function(targets,
						xcmsOptions,
                        libs="LipidPos",
                        RTs="none",
                        nCE=1,
                        corThresh=0.8,
                        checkIsotope=TRUE,
                        tolerance=25,
                        maxMZdiff=0.01,
                        matchWeight=0.5){
    
    ## Initialize results object-----------------------------------------------
    results <- initializeResultsAIF(targets)
    
    ## RTs intervals specified? --------------------------------------------
    if(RTs == "none") {
        message("No RT information provided...")
    } else if(is.data.frame(RTs)){
        message("Using user provided RT information...")
    } else stop("RTs must be a data.frame")
    
    ## load libraries ---------------------------------------------------------
    if(libs == "LipidPos") {
    	libraries <- MetaboAnnotatoR::LipidPos
    } else if(libs == "LipidNeg") {
    	libraries <- MetaboAnnotatoR::LipidNeg
    } else if(libs == "MetabolitesPos") {
    	libraries <- MetaboAnnotatoR::MetabolitesPos
    } else if(libs == "MetabolitesNeg") {
    	libraries <- MetaboAnnotatoR::MetabolitesNeg
    } else libraries <- loadLibs(libs)
    
    libfiles <- libraries$libfiles
    lib <- libraries$lib
    
    # evaluate if one or more samples are to be read---------------------------
    if(length(unique(targets[,3])) == 1){
        message("Reading and peak-picking data...")
        MSData <- readData(filePath=targets[1,3], 
                            xcmsOptions=xcmsOptions, nCE)
        xcmsF1 <- MSData$xcmsF1
        xcmsF2 <- MSData$xcmsF2
        peaksF1 <- MSData$peaksF1
        peaksF2 <- MSData$peaksF2
    } else NULL

    ## process each feature from the targets table-----------------------------
    for(i in seq_len(nrow(targets))){
        progNote <- paste("... Processing feature", i, "of", nrow(targets),
                            "...")
        message(progNote)
        if(length(unique(targets[,3]))>1){
            message("Reading and peak-picking data...")
            MSData <- readData(filePath=targets[i,3], 
                            xcmsOptions=xcmsOptions, nCE)
            xcmsF1 <- MSData$xcmsF1
            xcmsF2 <- MSData$xcmsF2
            peaksF1 <- MSData$peaksF1
            peaksF2 <- MSData$peaksF2
        } else NULL
    
        fmz <- targets[i,1]
        frt <- targets[i,2]
    
        ## get sample name to save files into --------------------------------
        SpName <- targets[i,3]
    
        ## get MS spectra at feature RT --------------------------------------
        message("Obtaining pseudo-MS/MS spectrum...")
        try(
            specs <- getPseudoMSMS(fmz, frt, xcmsF1, xcmsF2, peaksF1, peaksF2,
                            cthres1=corThresh, cthres2=corThresh)
        )
    
        if(exists("specs")){
            highCESpec <- specs$ms2
            pseudoSpec <- specs$aif
            inSourceSpec <- specs$insource
            ms2eic <- specs$ms2_eic
            feic <- specs$feic
            results$pseudoMSMS[[i]] <- specs
        } else results$pseudoMSMS[[i]] <- NA
    
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
        message("Matching candidate(s) fragments to pseudo-MS/MS spectra...")
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
                results$rankedResult[[i]] <- NA
                results$rankedSpectra[[i]] <- NA
            } else { rankedResult$pseudoMSMS <- "TRUE"
            }
        } else {
            output <- rankScore(result,specMatch)
            rankedResult <- output$rankedResult
            rankedSpec <- output$rankedSpecMatch
            results$rankedResult[[i]] <- rankedResult
            results$rankedSpectra[[i]] <- rankedSpec
        }

        ## Store highest rank annotation in global results---------------------
        if(!exists("rankedResult")) rankedResult <- NULL
        results$global[i,] <- storeAnnotations(global=results$global[i,], 
                                                rankedResult[1,])
    }
    # check ESI polarity
    polarity <- unique(polarity(xcmsF1))
    if(polarity == 1) polarity <- "positive"
    if(polarity == 0) polarity <- "negative"
    
    # store general options
    df <- data.frame(dataType="AIF", 
                        polarity=polarity, libraries=libs, nCE=nCE,
                        RTs=RTs, corThresh=corThresh,
                        checkIsotope=checkIsotope, matchWeight=matchWeight,
                        tolerance=paste(tolerance, "ppm"),
                        maxMZdiff=paste(maxMZdiff, "Da"), row.names="parameter")
    results$options <- as.data.frame(t(df))

    message("Job done!")
    return(results)
}

## Helper functions------------------------------------------------------------

## Initialization of results folder and global results table-------------------
initializeResultsAIF <- function(targets){
    # create objects to store results and metadata
    Date <- Sys.Date()
    Time <- format(Sys.time(), "%X")
    rankedResult <- list()
    rankedSpectra <- list()
    pseudoMSMS <- list()
    options <- NULL
    # create table to store global results
    global <- targets
    global[,c("metabolite", "feature.type", "ion.type", "isotope", 
                "mz.metabolite","matched.mz", "mz.error", "pseudoMSMS", 
                "fraction", "score")] <- NA
    # return global results table and results path as list
    results <- list(global=global,
    				Date=Date, 
    				Time=Time,
    				options=options,
    				rankedResult=rankedResult, 
    				rankedSpectra=rankedSpectra, 
    				pseudoMSMS=pseudoMSMS)
    return(results)
}

## read and peak-pick data from mzML or CDF files------------------------------
readData <- function(filePath, xcmsOptions, nCE){
    if(length(grep(".mzML", filePath)) == 1){
        # separate the two MS functions
        xcmsF1 <- MSnbase::readMSData(filePath, msLevel. = 1, mode = "onDisk")
        xcmsF2 <- MSnbase::readMSData(filePath, msLevel. = 2, mode = "onDisk")
        if(nCE > 1){
            maxCE <- max(xcmsF2@featureData@data$collisionEnergy)
            highCEscans <- which(
                xcmsF2@featureData@data$collisionEnergy == maxCE
                )
            xcmsF2@featureData@data <- xcmsF2@featureData@data[highCEscans,]
        }
    }
	# 2 CDF files have been used in Waters data, 
	# converted from raw using data bridge 
    if(length(grep(".CDF", filePath)) == 1){
        filePath2 <- grep("01.CDF", "02.CDF", filePath)
        # read the two MS functions
        xcmsF1 <- MSnbase::readMSData(filePath, mode = "onDisk")
        xcmsF2 <- MSnbase::readMSData(filePath2, mode = "onDisk")
    }
    # use same centwave parameters for both no- and high-collision scans
    cwp <- xcms::CentWaveParam(snthresh = xcmsOptions[3,2],
                                noise = xcmsOptions[2,2],
                                ppm = xcmsOptions[1,2],
                                peakwidth = as.numeric(xcmsOptions[4,2:3]),
                                prefilter = as.numeric(xcmsOptions[5,2:3]))
    peaksF1 <- xcms::findChromPeaks(xcmsF1, msLevel = 1L, param = cwp)
    peaksF2 <- xcms::findChromPeaks(xcmsF2, msLevel = 2L, param = cwp)
    # gather all data objects into a list
    MSData <- list(xcmsF1=xcmsF1, xcmsF2=xcmsF2, 
                    peaksF1=peaksF1, peaksF2=peaksF2)
    return(MSData)
}