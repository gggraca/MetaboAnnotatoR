#' Annotation of features using LC-MS AIF datasets processed using RAMClustR.
#'
#' Search possible matches of a feature in Lipid fragments and other
#' small molecule libraries.
#'
#' @author Goncalo Graca & Yuheng (Rene) Cai (Imperial College London)
#'
#' @param targetTable A csv file containing the list of features to annotate.
#' @param xcmsObject XCMS object containing the processed AIF datasets.
#' @param ramclustObj RAMClustR object with parent-fragment reconstructions
#' (clusters). See RAMClustR paper for more details
#' (https://pubs.acs.org/doi/10.1021/ac501530d).
#' @param ESImode Ionization mode: 'POS' for positive (default) or 'NEG'
#' for negative ionisation modes.
#' @param libs Fragment library to use: 'Lipids' (default) or 'Metabolites'
#' for other small molecules.
#' @param RTfile optional csv file with Lipid/metabolite classes Retention
#' Times in seconds
#' @param checkIsotope Whether or not to check the isotope type;
#' default is set to TRUE.
#' @param tolerance Tolerance in ppm for the candidate search.
#' @param maxMZdiff Maximum m/z difference between candidate fragments and
#' pseudo-MS/MS or AIF ions in Da.
#' @param matchWeight weight of the fragment matches to the final score;
#' value between 0 and 1; the remaining fraction of the weight comes from the
#' candidate m/z error.
#' @param ncandidates Maximum number of candidates to plot and store.
#' @return For each feature in the targeTable the will return a ranked list of
#' annotations, a plot pseudo-MS/MS spectrum for the matched ions,
#' a targeTable annotated with rank 1 annotations and a table with the options
#' used for the function.
#' @examples
#' # Run annotation of lipid features for positive LC-MS data XCMS set (XSet) 
#' # processed with RAMClustR:
#' annotateRC(targetTable="targetTable.csv", xcmsObject=XSet, ramclustObj=RC, 
#' libs="Lipids", ESImode="POS")
#' @importFrom utils read.csv write.csv
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics text
#' @export
annotateRC <- function(targetTable,
                       xcmsObject,
                       ramclustObj,
                       libs="Lipids",
                       ESImode="POS",
                       RTfile="none",
                       checkIsotope=TRUE,
                       tolerance=25,
                       maxMZdiff=0.01,
                       matchWeight=0.5,
                       ncandidates=5){
  
## Read targets table-------------------------------
targets <- read.csv(targetTable, header=TRUE)

## Initialize results object------------------------
results <- initializeResultsRC(targets, libs, ESImode)

## Get RTs intervals from file --------------------
RTs <- RTsFromFile(RTfile)

## load libraries ---------------------------------
libraries <- loadLibs(libs, ESImode)
libfiles <- libraries$libfiles
lib <- libraries$lib

## process each feature from the targets table------
for(i in seq_len(nrow(targets))){
    progNote <- paste("##### Processing feature", i, "of", nrow(targets), 
                      " ######")
    message(progNote)
    fmz <- targets[i,1]
    frt <- targets[i,2]

    # get MS spectra at feature RT --------------------
    #message("Reading high collision energy and pseudo-MS/MS spectra...")
    pseudoSpec <- RCspec(fmz, frt, ramclustObj)
    inSourceSpec <- xcmsSpec(fmz, frt, xcmsObject, highCE = FALSE)
    highCESpec <- xcmsSpec(fmz, frt, xcmsObject, highCE = TRUE)

    # Isotope check --------------------
    if(!checkIsotope) iso <- 0 else {
      #message("Checking isotope type...")
      iso <- checkIsotope(fmz, frt, inSourceSpec)
      }

	# Search Libraries --------------------
    if(is.null(pseudoSpec) & is.null(highCESpec)) { next
        } else {
    message("Searching candidates...")
    candidates <- searchLib(lib, libfiles, fmz - iso * 1.0034, frt,
	                        tolerance=tolerance, RTs, inSourceSpec)
    }

    # Compare fragments between Library candidates and high-collision-energy /
    # pseudo-MS/MS spectra --------------------
    if(is.null(pseudoSpec) & is.null(highCESpec) | length(unlist(candidates)) == 0) {
        result <- NULL
    } else {
        message("Matching candidate(s) fragments to pseudo-MS/MS and highCE spectra...")
        output <- mapply(compFrag, candidates,
                     lapply(as.numeric(names(candidates)),function(x) lib[[x]]),
                     MoreArgs=list(fmz, frt, iso, highCESpec, pseudoSpec,
                                     maxMZdiff=maxMZdiff,
                                     matchWeight=matchWeight), SIMPLIFY=FALSE)
        result <- do.call(rbind, lapply(output, "[[", 1))
        specMatch <- unlist(lapply(output, "[[", 2), recursive = FALSE)
        specMatch <- specMatch[!(specMatch) == "NULL"]
    }

    # Score ranking --------------------
    if(is.null(result)) {
      rankedResult <- targets[i,1:2]
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
        output <- rankScore(result, specMatch)
        rankedSpec <- output$rankedSpecMatch
        rankedResult <- output$rankedResult
    }
    
    ## Save output of ranked annotations -------
    if(exists("rankedResult")){
        saveRanked(fmz, frt, ncandidates, rankedResult, 
                   DatasetName=results$DatasetName, 
                   resultsDir=results$resultsDir)
    }

    ## Store highest rank annotation in global results------
    if(!exists("rankedResult")) rankedResult <- NULL
    results$global[i,] <- storeAnnotations(global=results$global[i,], 
                                          rankedResult[1,])
    
    ## save plot of matched spectra for the top n candidates
    if(exists("rankedSpec")){
        saveMatched(fmz, frt, highCESpec, result, output, ncandidates,
                    rankedSpec, DatasetName=results$DatasetName,
                    resultsDir=results$resultsDir)
    }
    ## Save pseudo-MS/MS spectrum per target feature
    if(exists("pseudoSpec")){
        saveRCpseudoMSMS(fmz, frt, pseudoSpec, 
                       DatasetName=results$DatasetName, 
                       resultsDir=results$resultsDir)
    }
}

## save global results table
write.csv(results$global, 
          file=paste(results$resultsDir, "Global_Results", ".csv", sep=""),
          row.names=FALSE)
## save general options
df <- data.frame(targetsTable_file=targetTable,
                 RAMClusterObject=as.character(substitute(ramclustObj)),
	             libraries=libs,
	             ESImode=ESImode,
	             RTfile=RTfile,
	             checkIsotope=checkIsotope,
	             matchWeight=matchWeight,
	             tolerance=paste(tolerance, "ppm"),
	             maxMZdiff=paste(maxMZdiff, "Da"),
	             row.names="Option")
df <- as.data.frame(t(df))
write.csv(df, file=paste(results$resultsDir, "General_options", ".csv", sep=""))
message('Job done!')

}


## Helper functions-------------------------------------------------------------

## Initialization of results folder and global results table
initializeResultsRC <- function(targets, libs, ESImode){
    # Create directory to store the results
    mainDir <- "./Annotations"
    Date <- Sys.Date()
    Time <- format(Sys.time(), "%X")
    Time <- gsub(":", "_", Time)
    subDir <- paste(libs, "_", ESImode, "_", "RAMClustR", "_",
                    Date, "_", Time, sep = "")
    dir.create(file.path(mainDir), showWarnings = FALSE)
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    resultsDir <- paste(mainDir, "/", subDir, "/", sep="")
    # create table to store global results
    global <- targets
    global[,c("metabolite", "feature.type", "ion.type", "isotope",
              "mz.metabolite", "matched.mz", "mz.error",
              "pseudoMSMS", "fraction", "score")] <- NA
    # create dataset name
    DatasetName <- paste(libs, "_", ESImode, "_", "RAMClustR", sep = "")
    # return global results table and results path as list
    results <- list(global=global, resultsDir=resultsDir, 
                    DatasetName=DatasetName)
    return(results)
}

## Save output of ranked annotations ------------------
saveRanked <- function(fmz, frt, ncandidates, rankedResult, 
                       DatasetName, resultsDir){
    write.csv(rankedResult[rankedResult$rank <= ncandidates,],
              file=paste(resultsDir, DatasetName, "_", round(fmz,3), "mz_", 
                           round(frt,3), "s_", "ranked_candidates.csv", 
                           sep=""), row.names=FALSE)
}

## Save plot of matched spectra for the top n candidates
saveMatched <- function(fmz, frt, 
                        highCESpec,
                        result,
                        output, 
                        ncandidates,
                        rankedSpec, 
                        DatasetName, 
                        resultsDir){
    if(!is.null(result)){
        if(length(rankedSpec) == 1) {
            plotCandidatesRC(fmz, frt, highCESpec, DatasetName, output, 1,
                             DirPath=resultsDir)
        }
        if(length(rankedSpec) <= ncandidates) {
            plots <- lapply(seq_len(length(rankedSpec)),
                            function(x) plotCandidatesRC(fmz,frt,
                                                         highCESpec,
                                                         DatasetName,
                                                         output, x,
                                                         resultsDir))
        }
        if(length(rankedSpec) > ncandidates) {
            plots <- lapply(seq_len(length(rankedSpec[1:ncandidates])),
                            function(x) plotCandidatesRC(fmz,frt,
                                                         highCESpec,
                                                         DatasetName,
                                                         output,x,
                                                         resultsDir))
        }
    }	else NULL
}

## Save pseudo-MS/MS spectrum per feature
saveRCpseudoMSMS <- function(fmz, frt, pseudoSpec, DatasetName, resultsDir){
    if(!is.null(pseudoSpec) & length(pseudoSpec) > 0){
        if (is.vector(pseudoSpec)) pseudoSpec <- as.data.frame(t(pseudoSpec))
        pdf(file = paste(resultsDir, "pseudoMSMS_",
                         DatasetName, "_", round(fmz,3), "mz_", round(frt,3),
                         "s", ".pdf", sep=""), width=8, height=5)
        
        plot(pseudoSpec[,1], pseudoSpec[,2], type='h',
             xlim=c(50, max(pseudoSpec[,1]) + 100),
             ylim=c(0, max(pseudoSpec[,2]) + max(pseudoSpec[,2])/1.5),
             xlab="m/z", ylab="intensity (a.u.)", col="black", lwd=1,
             main=paste("Pseudo-MS/MS Feature:", fmz, "m/z,", frt, "s"),
             cex.main=0.95, bty="L", xaxs="i", yaxs="i")
        text(pseudoSpec[,1] - 10, pseudoSpec[,2],
             as.character(round(pseudoSpec[,1],3)), pos=4, cex=0.8, srt=45)
        dev.off()
        write.csv(pseudoSpec, file=paste(resultsDir, "pseudoMSMS_",
                                         DatasetName, "_", round(fmz,3), "mz_",
                                         round(frt,3), "s", ".csv", sep = ""), 
                  row.names=FALSE)
    } else NULL
}
