#' @title Annotate features from LC-MS AIF datasets processed with 
#' RAMClustR
#'
#' @description
#' This function annotates features from LC-MS AIF using pseudo-MS/MS spectra 
#' obtained using RAMClustR, by matching ions to metabolite/lipid fragment
#' libraries.
#'
#' @author Goncalo Graca & Yuheng (Rene) Cai (Imperial College London)
#'
#' @param targets A data.frame containing the features to annotate
#' and the file paths to the raw data.
#' @param xcmsObject XCMS object containing the processed AIF datasets.
#' @param ramclustObj RAMClustR object with parent-fragment reconstructions
#' (clusters). See RAMClustR paper for more details
#' (https://pubs.acs.org/doi/10.1021/ac501530d).
#' @param libs Fragment libraries to use. Either the built-in libraries can be
#' specified (\code{LipidPos}, \code{LipidNeg}, \code{MetabolitesPos}, 
#' \code{MetabolitesNeg}) or the full path to user-defined libraries.
#' @param RTs Optional data.frame with Lipid/metabolites classes Retention
#' Times in seconds.
#' @param checkIsotope Whether or not to check the isotope type;
#' default is set to TRUE.
#' @param tolerance Tolerance in ppm for the candidate search.
#' @param maxMZdiff Maximum m/z difference between candidate fragments and
#' pseudo-MS/MS or AIF ions in Da.
#' @param matchWeight weight of the fragment matches to the final score;
#' value between 0 and 1; the remaining fraction of the weight comes from the
#' candidate m/z error.
#' @return For each feature in the targets data frame the function will return 
#' a list containing: a data frame with with rank 1 annotations 
#' (\code{global}), the the date and time of annotation, a data frame with the 
#' annotation options. For each feature the following lists are returned: 
#' ranked annotations for each feature (\code{rankedResults}), 
#' the corresponding ranked matched spectra (\code{rankedSpectra}), 
#' the pseudo-MS/MS spectra (\code{pseudoMSMS}), in-source spectra 
#' (\code{inSourceSpectra}) and AIF spectrum (\code{AIFspectra}).
#' @examples
#' # Read RAMClustR (RC) and XCMS processed example data lipid positive 
#' # LC-MS data
#' data("RC")
#' data("xset")
#' # read the table containing features to annotate
#' tfile <- system.file("extdata", "targetTable.csv", 
#' package="MetaboAnnotatoR")
#' targets <- read.csv(tfile)
#' # Run the annotation procedure
#' annotations <- annotateRC(targets, xcmsObject=xset, ramclustObj=RC, 
#' libs="LipidPos", RTs="none", checkIsotope=TRUE)
#' @export
annotateRC <- function(targets, xcmsObject, ramclustObj,
                        libs="LipidPos", RTs="none",
                        checkIsotope=TRUE, tolerance=25,
                        maxMZdiff=0.01, matchWeight=0.5){

    ## Initialize results object------------------------
    results <- initializeResultsRC(targets)

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

    ## process each feature from the targets table------
    for(i in seq_len(nrow(targets))){
        progNote <- paste("... Processing feature", i, "of", nrow(targets), 
                            "...")
        message(progNote)
        fmz <- targets[i,1]
        frt <- targets[i,2]

        # get MS spectra at feature RT --------------------
        #message("Reading high collision energy and pseudo-MS/MS spectra...")
        pseudoSpec <- RCspec(fmz, frt, ramclustObj)
        inSourceSpec <- xcmsSpec(fmz, frt, xcmsObject, highCE = FALSE)
        highCESpec <- xcmsSpec(fmz, frt, xcmsObject, highCE = TRUE)
        # store spectra in the annotation results object
        if(!is.null(pseudoSpec)) {
            results$pseudoMSMS[[i]] <- pseudoSpec
            } else results$pseudoMSMS[[i]] <- NA
        results$AIFspectra[[i]] <- highCESpec
        results$inSourceSpectra[[i]] <- inSourceSpec

        # Isotope check --------------------
        if(!checkIsotope) iso <- 0 else {
            iso <- checkIsotope(fmz, frt, inSourceSpec)
        }

        # Search Libraries --------------------
        if(is.null(pseudoSpec) & is.null(highCESpec)) { next
        } else {
            message("Searching candidates...")
            candidates <- searchLib(lib, libfiles, fmz - iso * 1.0034, frt,
                                    tolerance=tolerance, RTs, inSourceSpec)
        }
        #Compare fragments between Library candidates and high-collision-energy
        #pseudo-MS/MS spectra --------------------
        if(is.null(pseudoSpec) & 
            is.null(highCESpec) | 
            length(unlist(candidates)) == 0){
            result <- NULL
        } else {
            message("Matching fragments to pseudo-MS/MS and highCE spectra...")
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

            # type of ion isotope
            rankedResult$isotope <- paste("M+", iso, sep="")

            # pseudoMSMS flag
            if(is.null(pseudoSpec)) {
                rankedResult$pseudoMSMS <- "FALSE"
            } else rankedResult$pseudoMSMS <- "TRUE"

        } else {
            output <- rankScore(result, specMatch)
            rankedSpec <- output$rankedSpecMatch
            rankedResult <- output$rankedResult
            results$rankedResult[[i]] <- rankedResult
            results$rankedSpectra[[i]] <- rankedSpec
        }
        ## Store highest rank annotation in global results------
        results$global[i,] <- storeAnnotations(global=results$global[i,], 
                                rankedResult[1,])
    }
    # store options
    df <- data.frame(dataType="RAMClustR",
                    polarity=ramclustObj$ExpDes$instrument$value[9],
                    libraries=libs,
                    RTs=RTs,
                    checkIsotope=checkIsotope, matchWeight=matchWeight,
                    tolerance=paste(tolerance, "ppm"),
                    maxMZdiff=paste(maxMZdiff, "Da"), row.names="parameter")
    results$options <- as.data.frame(t(df))

    message('Job done!')

    return(results)
}


## Helper functions------------------------------------------------------------

## Initialization of results folder and global results table
initializeResultsRC <- function(targets){
    # create objects to store results and metadata
    Date <- Sys.Date()
    Time <- format(Sys.time(), "%X")
    rankedResult <- list()
    rankedSpectra <- list()
    pseudoMSMS <- list()
    AIFspectra <- list()
    inSourceSpectra <- list()
    options <- NULL
    # create table to store global results
    global <- targets[,c(1,2)]
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
                    pseudoMSMS=pseudoMSMS,
                    AIFspectra=AIFspectra,
                    inSourceSpectra=inSourceSpectra)
    return(results)
}
