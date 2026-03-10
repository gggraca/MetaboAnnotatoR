#' @title Save annotation results
#'
#' @description
#' Saves all annotation related data, such as annotation table, options,
#' and optionally plots of pseudo-MS/MS and matched fragments spectra. 
#'
#' @author Goncalo Graca & Yuheng (Rene) Cai (Imperial College London)
#'
#' @param annotations Annotation object created from running annotation
#' \code{annotateAIF} or \code{annotateRC} functions.
#' @param saveOptions If \code{TRUE}, will save the annotation options as 
#' .csv file.
#' @param saveXCMSoptions Saves the XCMS options if the annotations
#' originate from AIF raw files, if set to \code{TRUE}.
#' @param saveRanked Option to save the ranked candidate table and 
#' respective pseudo-MS/MS. The default option is \code{TRUE}.
#' @param savePseudoMSMS Option to save the pseudo-MS/MS for the features 
#' from the targets table. The default option is code{FALSE}. If \code{TRUE},
#' the pseudo-MS/MS will be saved as \code{.png} images and \code{.mgf} files.
#' Note that EICs will only be saved in the png image if the annotations
#' originate from AIF raw files.
#' @param DirPath Path to the folder where the plots (as \code{.png}), 
#' \code{.csv} tables and \code{.mgf} files will be saved.
#' @return Global and candidate annotations as .csv files and pseudo-MS/MS 
#' spectra as \code{.png} and/or \code{.mgf} files.
#' @examples
#' #' # Read RAMClustR (RC) and XCMS processed example data:
#' tfile <- system.file("targetTable.csv", package = "MetaboAnnotatoR")
#' targets <- read.csv(tfile)
#' # Run annotation of lipid features for positive LC-MS 
#' # processed with RAMClustR:
#' annotations <- annotateRC(targets, xcmsObject=xset, ramclustObj=RC, 
#' libs="LipidPos", RTs="none", checkIsotope=TRUE)
#' # Finally, save the results to a user-defined directory
#' userDir <- tempdir()
#' saveAnnotations(annotations, DirPath=userDir, saveOptions=TRUE, 
#' saveXCMSoptions=FALSE, saveRanked=TRUE,
#' saveRankedSpec=TRUE, savePseudoMSMS=TRUE)
#' @importFrom utils read.csv write.csv
#' @export
saveAnnotations <- function(annotations,
							DirPath="",
							saveOptions=TRUE, 
							saveXCMSoptions=FALSE, 
							saveRanked=TRUE,
							saveRankedSpec=FALSE,
							savePseudoMSMS=FALSE){
	
	# separate useful variables from the annotations list
	global <- annotations$global
	options <- annotations$options
	dataType <- options[1,]
	polarity <- options[2,]
	rankedResult <- annotations$rankedResult
	rankedSpectra <- annotations$rankedSpectra
	pseudoMSMS <- annotations$pseudoMSMS
	
	# convert polarity label for saveMgf
	if(polarity == "positive") polarity <- as.integer(1)
	if(polarity == "negative") polarity <- as.integer(0) 
	
	# save options
	if(saveOptions){
		fpath <- file.path(DirPath, "annotation_options.csv")
		write.csv(options, fpath)
	}
	
	# save global annotations result to file
	fpath <- file.path(DirPath, "annotation_global_results.csv")
	write.csv(global, fpath, row.names=FALSE)
		
	# save ranked results
	if(saveRanked){
		for(i in seq_along(rankedResult)){
			if(!is.null(rankedResult[[i]])){
				fpath <- file.path(DirPath, "/", 
								   round(global$feature.mz[i], 3),
								   "mz_",
								   round(global$feature.rt[i], 3),
								   "s_", 
								   "ranked_results.csv",
								   fsep="")
				write.csv(rankedResult[[i]], fpath, , row.names=FALSE)
			}
		}
	}

    # save ranked spectra
	if(saveRankedSpec){
		if(dataType == "RAMClustR"){
			l <- nrow(global)
			l <- seq_len(l)
			l <- lapply(l, 
						function(x) plotCandidatesRC(
							rankedResult[[x]], 
							rankedSpectra[[x]], 
							DirPath)
						)
		}
		if(dataType == "AIF") {
			l <- nrow(global)
			l <- seq_len(l)
			l <- lapply(l, 
						function(x) plotCandidatesAIF(
							rankedResult[[x]], 
							rankedSpectra[[x]],
							pseudoMSMS[[x]],
							DirPath)
			)
		}
	}
	
	# save pseudoMSMS as .mgf
	if(savePseudoMSMS & dataType == "AIF"){
		saveAIFmgf(global, pseudoMSMS, polarity, DirPath)
	}
	if(savePseudoMSMS & dataType == "RC"){
		saveRCmgf(global, pseudoMSMS, polarity, DirPath)
	}
}

## helper functions------------------------------------------------------------
# save pseudo-MS/MS spectrum as .mgf format------------------------------------
saveRCmgf <- function(global, pseudoMSMS, polarity, DirPath){
	for(i in seq_along(pseudoMSMS)){
		if(!is.na(pseudoMSMS)){
			spec <- new("Spectrum2", 
						mz=pseudoMSMS[[i]][,"mz"], 
						intensity=pseudoMSMS[[i]][,"into"],
						precursorMz=global$feature.mz[i],
						rt=global$feature.rt[i],
						polarity=polarity,
						centroided=TRUE)
			fpath <- file.path(DirPath, "/", 
							   round(global$feature.mz[i], 3),
							   "mz_",
							   round(global$feature.rt[i], 3),
							   "s_", 
							   "pseudo-MSMS", ".mgf",
							   fsep="")
			MSnbase::writeMgfData(spec, fpath)
		} else next
	}
}


saveAIFmgf <- function(global, pseudoMSMS, polarity, DirPath){
	for(i in seq_along(pseudoMSMS)){
		if(!is.null(pseudoMSMS[[i]]$aif)){
			spec <- new("Spectrum2", 
						mz=pseudoMSMS[[i]]$aif[,"mz"], 
						intensity=pseudoMSMS[[i]]$aif[,"into"],
						precursorMz=global$feature.mz[i],
						rt=global$feature.rt[i],
						polarity=polarity,
						centroided=TRUE)
			fpath <- file.path(DirPath, "/", 
							   round(global$feature.mz[i], 3),
							   "mz_",
							   round(global$feature.rt[i], 3),
							   "s_", 
							   "pseudo-MSMS", ".mgf",
							   fsep="")
			MSnbase::writeMgfData(spec, fpath)
		} else next
	}
}

# plotCandidatesRC to plot annotation candidates matched spectra--------------
plotCandidatesRC <- function(rankedResult, rankedSpectra, DirPath){
	if(!is.null(rankedResult)){
		candidates <- nrow(rankedResult)
		if(!is.null(candidates)) candidates <- seq_len(candidates)
		
		l <- lapply(candidates, function(x) plotSingleCandidateRC(rankedResult, 
														rankedSpectra, 
														pseudoMSMS,
														x,
														DirPath))
	}
}

# plot one RC candidate and save to DirPath
plotSingleCandidateRC <- function(rankedResult, rankedSpectra, pseudoMSMS,
								  candidate, DirPath){
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
	
	fpath <- file.path(round(fmz,3), "mz_", 
					   round(frt,3), "s", "_candidate_", candidate,".pdf", fsep="")
	ggplot2::ggsave(fpath, plot=plt, path=DirPath, width=10, height=5)
}

# plotCandidatesAIF to plot annotation candidates matched spectra--------------
plotCandidatesAIF <- function(rankedResult, 
							  rankedSpectra, 
							  pseudoMSMS,
							  DirPath){
	candidates <- nrow(rankedResult)
	if(!is.null(candidates)) candidates <- seq_len(candidates)
	
	l <- lapply(candidates, function(x) plotSingleCandidateAIF(rankedResult, 
                                                 rankedSpectra, 
												 pseudoMSMS,
												 x,
												 DirPath))
}

# plot one AIF candidate and save to DirPath
plotSingleCandidateAIF <- function(rankedResult, 
								   rankedSpectra, 
								   pseudoMSMS,
								   candidate,
								   DirPath){
	# get relevant information from the rankedCResult and rankedSpectra objects
	metabolite <- rankedResult[candidate,"metabolite"]
	ionType <- rankedResult[candidate,"ion.type"]
	score <- rankedResult[candidate,"score"]
	adductName <- paste(metabolite,ionType)
	candidateMZ <- round(rankedResult[candidate, "mz.metabolite"],3)
	MZerror <- round(rankedResult[candidate,"mz.error"],1)
	specMatch <- rankedSpectra[[candidate]]
	highCESpec <- pseudoMSMS$ms2
	ms2eic <- pseudoMSMS$ms2_eic
	rnk <- rankedResult[candidate,"rank"]
	fmz <- rankedResult[candidate,"feature.mz"]
	frt <- rankedResult[candidate,"feature.rt"]
	
	if(nrow(specMatch) >= 1){
		# plotting part
		# EICs
		mz_idx <- match(specMatch[,1], highCESpec[,"mz"])
		df1 <- lapply(mz_idx, function(x) data.frame(
			intensity=intensity(ms2eic[x,1]), rt=rtime(ms2eic[x,1]),
			mz=as.character(rep(paste(round(highCESpec[x,"mz"],3),"m/z")))))
		df1 <- do.call("rbind", df1)
		p1 <- ggplot2::ggplot(df1[!is.na(df1$intensity),],
							  ggplot2::aes(x=rt, y=intensity, colour=mz)) +
			ggplot2::geom_point() +
			ggplot2::geom_line() +
			ggplot2::labs(x="RT (s)", y = "Intensity (a.u.)", 
						  colour="fragments") +
			ggplot2::ggtitle(paste("Feature",round(fmz,3),"m/z",round(frt),
								   "s,","Rank", rnk, "result:", adductName,
								   ", mz.error =", MZerror, "ppm", ", 
                                    score =", round(score, 2))) +
			ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
		# Spectrum
		df2 <- as.data.frame(specMatch[,c("mz", "into")])
		p2 <- ggplot2::ggplot(df2,
							  ggplot2::aes(x=mz, y=into, label=round(mz, 3))) +
			ggplot2::geom_segment(ggplot2::aes(xend=mz, yend=0),
								  color="red", lwd=0.5) +
			ggplot2::geom_text(size=3, angle=45, hjust=0, vjust=0) +
			ggplot2::ggtitle(paste("ions matched to", adductName)) +
			ggplot2::theme_minimal() +
			ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
			ggplot2::ylim(0, max(df2[,2]) + 0.1*max(df2[,2])) +
			ggplot2::xlim(min(df2[,1])-50, max(df2[,1])+50) +
			ggplot2::labs(x = "m/z", y = "Intensity (a.u.)")
		fpath <- file.path(round(fmz,3),"mz_", round(frt,3), "s", "_candidate_",
						   candidate, ".pdf", fsep="")
		plt <- gridExtra::grid.arrange(p1, p2, nrow=2)
		ggplot2::ggsave(fpath, plot=plt, device="pdf", 
						path=DirPath, width=10, height=10)
	} else NULL
}