## Common helper functions to organize annotation results 
## Goncalo Graca (Imperial College London)

## Load libraries and get libraries filepaths----------------------------------
loadLibs <- function(libs){
	message("Loading user-defined libraries...")
	libfiles <- list.files(path=libs, full.names=TRUE)
	# check.names=FALSE to use the original header names in the annotations:
	lib <- lapply(libfiles, read.csv, header=TRUE, sep=",",
				  check.names=FALSE)
	libraries <- list(lib=lib, libfiles=libfiles)
	return(libraries)
}

## Store result for high rank candidate on global results table----------------
storeAnnotations <- function(global, rankedResult){
    if(!is.null(rankedResult)){
        global$isotope <- as.character(rankedResult$isotope)
        global$metabolite <- as.character(rankedResult$metabolite)
        global$mz.metabolite <- rankedResult$mz.metabolite
        global$matched.mz <- rankedResult$matched.mz
        global$mz.error <- rankedResult$mz.error
        global$ion.type <- as.character(rankedResult$ion.type)
        global$feature.type <- as.character(rankedResult$feature.type)
        global$pseudoMSMS <- rankedResult$pseudoMSMS
        global$fraction <- as.character(rankedResult$fraction)
        global$score <- rankedResult$score
    } else {
        global[,c("metabolite", "feature.type", "ion.type", "isotope",
                    "mz.metabolite", "matched.mz", "mz.error", "pseudoMSMS",
                    "fraction", "score")] <- NA
    }
    return(global)
}