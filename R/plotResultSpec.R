#' Plot pseudo-MS/MS composed of candidate matched ions
#' 
#' Function to visualise the spectra containing the matched ions to each
#' candidate annotation result
#' 
#' @author Goncalo Graca (Imperial College London)
#'
#' @param annotations An output object from the \code{annotateAIF} or 
#' \code{annotaterRC} functions.
#' @param feature Index of the annotated feature, as specified in the
#' \code{targets} data frame.
#' @param candidate Index of the candidate annotation as specified in the
#' \code{annotations$rankedResults} data frame.
#' @return A pseudo-MS/MS spectrum is plotted.
#' @examples
#' # Read RAMClustR (RC) and XCMS processed example data:
#' tfile <- system.file("targetTable.csv", package = "MetaboAnnotatoR")
#' targets <- read.csv(tfile)
#' # Run annotation of lipid features for positive LC-MS 
#' # processed with RAMClustR:
#' annotations <- annotateRC(targets, xcmsObject=xset, ramclustObj=RC, 
#' libs="LipidPos", RTs="none", checkIsotope=TRUE)
#' # plot the rank 1 candidate of the 3rd feature in the annotation$targets
#' plotResultSpec(annotations, 3, 1)
#' @export
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
