## Common helper functions to organize annotation results 
## Goncalo Graca (Imperial College London)
## These functions are not exported

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

## Plot pseudo-MS/MS composed of candidate matched ions
plotResultSpec <- function(annotations, feature, candidate){
	    rankedResult <- annotations$rankedResult[[feature]]
	    if(is.null(rankedResult)) {
	    	message("No annotation result found for this feature.")
	    } else {
	    	rankedSpectra <- annotations$rankedSpectra[[feature]]
	    	# get relevant information from the rankedResult and rankedSpectra objects
	    	metabolite <- rankedResult[candidate,"metabolite"]
	    	ionType <- rankedResult[candidate,"ion.type"]
	    	score <- rankedResult[candidate,"score"]
	    	adductName <- paste(metabolite, ionType)
	    	candidateMZ <- round(rankedResult[candidate,"mz.metabolite"], 3)
	    	MZerror <- round(rankedResult[candidate,"mz.error"], 1)
	    	rnk <- rankedResult[candidate,"rank"]
	    	specMatch <- rankedSpectra[[candidate]]
	    	fmz <- rankedResult[candidate,"feature.mz"]
	    	frt <- rankedResult[candidate,"feature.rt"]
	    	
	    	# plotting part
	    	if(is.null(nrow(specMatch))) {
	    		df <- as.data.frame(specMatch[c("mz", "into")])
	    	} else df <- as.data.frame(specMatch[,c("mz", "into")])
	    	
	    	plt <- ggplot2::ggplot(df,
	    						   ggplot2::aes(x=mz, y=into, label=round(mz, 3))) +
	    		ggplot2::geom_segment(ggplot2::aes(xend=mz, yend=0), 
	    							  color="red", lwd=0.5) +
	    		ggplot2::geom_text(size=3, angle=45, hjust=0, vjust=0) +
	    		ggplot2::ggtitle(paste("Feature ", round(fmz, 4), "m/z", "_", 
	    							   round(frt, 4), "s ", ", Rank ", rnk, 
	    							   " result: ", adductName, ", mz.error = ", 
	    							   MZerror, " ppm, ", "score = ", 
	    							   round(score, 2), sep="")) +
	    		ggplot2::theme_minimal() +
	    		ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +
	    		ggplot2::ylim(0, max(df$into) + 0.1*max(df$into)) +
	    		ggplot2::xlim(min(df$mz)-50, max(df$mz)+50) +
	    		ggplot2::labs(x="m/z", y="Intensity (a.u.)")
	    	plt
	    }
}
