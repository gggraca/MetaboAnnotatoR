#' @title Re-annotate features (MS/MS or pseudo-MS/MS)
#'
#' @description
#' This function re-annotates pseudo-MS/MS spectra from a MetaboAnnotatoR
#' annotation object (obtained from AIF data or ISF(MS1)) 
#' stored in the Global Environment. The spectra are generally acquired from 
#' an LC-MS experiment. Re-annotation is most useful when repeating annotation
#' using different fragment libraries or different candidate search and/or
#' scoring parameters.
#'
#' @author Goncalo Graca (Imperial College London)
#'
#' @param AnnotationObj A MetaboAnnotatoR object, containing the results from 
#' a previous annotation job.
#' @param libs Fragment libraries to use. Specify one of default libraries
#' provided as data object with the package (\code{LipidPos}, \code{LipidNeg}, 
#' \code{MetabolitesPos} and \code{MetabolitesNeg}) or the full path to 
#' user-defined libraries.
#' @param RTs Optional data.frame with Lipid/metabolites classes Retention
#' Times in seconds.
#' @param checkIsotope Whether or not to check the isotope type;
#' default is set to FALSE as the precursor ions will likely be the 
#' monoisotopic mass of the isotopic series.
#' @param tolerance Tolerance in ppm for the candidate search.
#' @param maxMZdiff Maximum m/z difference between candidate fragments and
#' pseudo-MS/MS or AIF ions in Da.
#' @param matchWeight weight of the fragment matches to the final score;
#' value between 0 and 1; the remaining fraction of the weight comes from the
#' candidate m/z error.
#' @return For each feature in the annotation object (\code{AnnotationObj}) 
#' the function will return a list containing: a data frame with with rank 1 
#' annotations (\code{global}), the the date and time of annotation, a 
#' data frame with the annotation options. For each feature the following 
#' lists are returned: ranked annotations for each feature 
#' (\code{rankedResults}), the corresponding ranked matched spectra 
#' (\code{rankedSpectra}), the pseudo-MS/MS spectra (\code{pseudoMSMS}), 
#' in-source spectra (\code{inSourceSpectra}) and AIF spectrum 
#' (\code{AIFspectra}, if available).
#' @examples
#' # Set a directory to save the example .mzML file
#' userDir <- tempdir()
#' # Download the example .mzML from zenodo website into the specified 
#' #directory as "Lipid_Positive_QC.mzML"...define file.path first:
#' fpath <- file.path(userDir, "Lipid_Positive_QC.mzML")
#' download.file(
#' "https://zenodo.org/records/17408169/files/Lipid_Positive_QC.mzML?download=1", 
#' fpath)
#' # create a new targetTable with one feature to annotate
#' # the Sample.name is the path to the mzML file
#' targets <- data.frame(feature.mz=520.3408533, feature.rt=100.6238759, 
#' Sample.name=fpath)
#' # read the default xcms parameters on the XCMS_options.csv file and modify
#' # the noise threshold parameter
#' xcmsOptionsPath <- system.file("extdata", "XCMS_options.csv", 
#' package="MetaboAnnotatoR")
#' xcmsOptions <- read.csv(xcmsOptionsPath)
#' xcmsOptions[2,2] <- 1000
#' # Read the default lipid positive libraries
#' data("LipidPos")
#' # Run the annotation using the built-in lipid POS library:
#' annotations <- annotateISF(targets, xcmsOptions, 
#' libs="LipidPos", RTs="none", corThresh=0.8,
#' checkIsotope=TRUE)
#' # Read the default metabolites positive libraries
#' data("MetabolitesPos")
#' # re-annotate using the metabolites positive libraries
#' newAnnotations <- reannotate(annotations, libs="MetabolitesPos", 
#' RTs="none", checkIsotope=FALSE)
#' @export
reannotate <- function(AnnotationObj,
                            libs="LipidPos", RTs="none",
                            checkIsotope=TRUE, tolerance=25,
                            maxMZdiff=0.01, matchWeight=0.5){
	
    ## organise spectra and related data from the AnnotationObj----------------
	# targets list
	targets <- AnnotationObj$global[,1:3]
    
    ## Initialize results object-----------------------------------------------
    results <- reinitializeResults(AnnotationObj)
    
    ## identify data type------------------------------------------------------
    dataType <- AnnotationObj$options[1,1]
    
    ## RTs intervals specified? -----------------------------------------------
    if(RTs == "none") {
        message("No RT information provided...")
    } else if(is.data.frame(RTs)){
        message("Using user provided RT information...")
    } else stop("RTs must be a data.frame")

    ## load libraries ---------------------------------------------------------
    if(libs == "LipidPos") {
        if(!exists("LipidPos")) {
            stop("LipidPos not found, please use data(LipidPos)")
        } else libraries <- LipidPos
    } else if(libs == "LipidNeg") {
        if(!exists("LipidNeg")) {
            stop("LipidNeg not found, please use data(LipidNeg)")
        } else libraries <- LipidNeg
    } else if(libs == "MetabolitesPos") {
        if(!exists("MetabolitesPos")) {
            stop("MetabolitesPos not found, please use data(MetabolitesPos)")
        } else libraries <- MetabolitesPos
    } else if(libs == "MetabolitesNeg") {
        if(!exists("MetabolitesNeg")) {
            stop("MetabolitesNeg not found, please use data(MetabolitesNeg)")
        } else libraries <- MetabolitesNeg
    } else libraries <- loadLibs(libs)
    
    libfiles <- libraries$libfiles
    lib <- libraries$lib

    ## process each feature from the targets table-----------------------------
    for(i in seq_len(nrow(targets))){
        progNote <- paste("... Processing feature", i, "of", nrow(targets), 
                            "...")
        message(progNote)
        fmz <- targets[i,1]
        frt <- targets[i,2]

        # read MS/MS spectra --------------------------------------------------
        if(dataType == "AIF"){
        	inSourceSpec <- results$pseudoMSMS[[i]]$insource
            pseudoSpec <- results$pseudoMSMS[[i]]$aif
            highCESpec <- results$pseudoMSMS[[i]]$ms2
        } else if(dataType == "ISF"){
        	inSourceSpec <- results$pseudoMSMS[[i]]$insource
            pseudoSpec <- results$pseudoMSMS[[i]]$insource
            highCESpec <- results$pseudoMSMS[[i]]$ms1
        }
        
        # Isotope check -------------------------------------------------------
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

        # Search Libraries ----------------------------------------------------
        if(is.null(pseudoSpec)) { next
        } else {
            message("Searching candidates...")
            candidates <- searchLib(lib, libfiles, fmz-iso*1.0034, frt,
                                    tolerance=tolerance, RTs, inSourceSpec)
        }
        #Compare fragments between Library candidates and high-collision-energy
        #pseudo-MS/MS spectra --------------------
        if(is.null(pseudoSpec) & 
            is.null(highCESpec) | 
            length(unlist(candidates)) == 0){
            result <- NULL
        } else {
            message("Matching fragments to pseudo-MS/MS...")
            output <- mapply(compFrag, candidates, 
                                lapply(as.numeric(names(candidates)),
                                function(x) lib[[x]]),
                                MoreArgs=list(fmz, frt, iso, highCESpec, 
                                                pseudoSpec,
                                                maxMZdiff=maxMZdiff,
                                                matchWeight=matchWeight), 
                                SIMPLIFY=FALSE)
            result <- do.call(rbind, lapply(output, "[[", 1))
            specMatch <- unlist(lapply(output, "[[", 2), recursive = FALSE)
            specMatch <- specMatch[!(specMatch) == "NULL"]
        }
        # Score ranking --------------------
        if(is.null(result)) {
            rankedResult <- targets[i,c(1,2)]
            rankedResult[c("metabolite", "feature.type", "ion.type", "isotope",
                            "mz.metabolite", "matched.mz", "mz.error", 
                            "pseudoMSMS", "fraction", "score")] <- NA
            rankedSpec <- NULL

        } else {
            output <- rankScore(result, specMatch)
            rankedSpec <- output$rankedSpecMatch
            rankedResult <- output$rankedResult
            # type of ion isotope
            rankedResult$isotope <- paste("M+", iso, sep="")
            # pseudoMSMS flag
            rankedResult$pseudoMSMS <- "TRUE"
            results$rankedResult[[i]] <- rankedResult
            results$rankedSpectra[[i]] <- rankedSpec
        }
        ## Store highest rank annotation in global results------
        results$global[i,] <- storeAnnotations(global=results$global[i,], 
                                rankedResult[1,])
    }
    # store options
    df <- data.frame(dataType=dataType,
                    polarity=AnnotationObj$options[2,1],
                    libraries=libs,
                    RTs=RTs,
                    checkIsotope=checkIsotope, matchWeight=matchWeight,
                    tolerance=paste(tolerance, "ppm"),
                    maxMZdiff=paste(maxMZdiff, "Da"),
    	            reannotation=TRUE, row.names="parameter")
    results$options <- as.data.frame(t(df))

    message('Job done!')

    return(results)
}


## Helper functions------------------------------------------------------------

## Initialization of results folder and global results table
reinitializeResults <- function(AnnotationObj){
    # create objects to store results and metadata
    Date <- Sys.Date()
    Time <- format(Sys.time(), "%X")
    rankedResult <- list()
    rankedSpectra <- list()
    pseudoMSMS <- AnnotationObj$pseudoMSMS
    options <- NULL
    # create table to store global results
    global <- AnnotationObj$global[,1:3]
    global[,c("metabolite", "feature.type", "ion.type", "isotope",
                "mz.metabolite", "matched.mz", "mz.error",
                "pseudoMSMS", "fraction", "score")] <- NA
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
