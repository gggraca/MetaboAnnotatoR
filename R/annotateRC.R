#' Annotation of features using LC-MS AIF datasets processed using RAMClusteR.
#'
#' Search possible matches of a feature in Lipid fragments and other
#' small molecule libraries.
#'
#' @author Goncalo Graca & Yuheng (Rene) Cai (Imperial College London)
#'
#' @param targetTable A csv file containing the list of features to annotate.
#' @param xcmsObject XCMS object containing the processed AIF datasets.
#' @param ramclustObj RAMClustR object with parent-fragment reconstructions.
#' (clusters).
#' @param ESImode Ionization mode: 'POS' for positive (default) or 'NEG'
#' for negative ionisation modes.
#' @param libs Fragment library to use: 'Lipids' (default) or 'Metabolites'
#' for other small molecules.
#' @param RTfile optional csv file with Lipid/metabolites classes Retention
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
#' @export
annotateRC <- function(targetTable,
                       xcmsObject,
                       ramclustObj,
                       libs = "Lipids",
                       ESImode = "POS",
                       RTfile = "none",
                       checkIsotope = TRUE,
                       tolerance = 25,
                       maxMZdiff = 0.01,
                       matchWeight = 0.5,
                       ncandidates = 5){
# Read targets table
targets <- read.csv(targetTable, header=TRUE)

# create table to store global results
global <- targets
global[,c("metabolite", "feature.type", "ion.type", "isotope",
          "mz.metabolite", "matched.mz", "mz.error",
          "pseudoMSMS", "fraction", "score")] <- NA

# Get RTs from file --------------------
	if(RTfile=="none") {
		message("No RT information provided...")
		RTs <- "none"
		} else {
		message("Reading RT information...")
		RTs <- read.csv(RTfile,header=TRUE)
	}

# Create directory to store the results
mainDir <- "./Annotations"
Date <- Sys.Date()
Time <- format(Sys.time(), "%X")
Time <- gsub(":", "_", Time)
subDir <- paste(libs, "_", ESImode, "_", "RAMClustR", "_",
                Date, "_", Time, sep = "")
dir.create(file.path(mainDir), showWarnings = FALSE)
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

# load libraries --------------------
if(file.exists("./Libraries/")){
  message("Loading user-defined libraries...")
  libfiles <- list.files(path = paste("./Libraries/",libs,"/",
                                      ESImode, sep = ""), full.names = TRUE)
} else {
  message("Loading default libraries...")
  defaultLibPath <- system.file(paste("/Libraries/",libs,"/",ESImode, sep = ""),
                                package = "MetaboAnnotatoR")
  libfiles <- list.files(defaultLibPath, full.names = TRUE)
}

# check.names=FALSE to use the original header names in the annotations:
lib <- lapply(libfiles, read.csv, header = TRUE, sep=",", check.names = FALSE)


for (i in 1:dim(targets)[1]){
	message(paste("####### Processing feature", i, "of",
	              dim(targets)[1])," ########")

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
	                        tolerance = tolerance, RTs, inSourceSpec)
	}

	# Compare fragments between Library candidates and high-collision-energy /
	# pseudo-MS/MS spectra --------------------
	if(is.null(pseudoSpec) & is.null(highCESpec) |
	   length(unlist(candidates))== 0) {
	  result <- NULL
	} else {
	message("Mathing candidate(s) fragments to pseudo-MS/MS and highCE spectra...")
	output <- mapply(compFrag, candidates,
	                 lapply(as.numeric(names(candidates)),function(x) lib[[x]]),
	                 MoreArgs = list(fmz, frt, iso, highCESpec, pseudoSpec,
                                 maxMZdiff = maxMZdiff,
                                 matchWeight = matchWeight),
	                 SIMPLIFY = FALSE)

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
		if (iso == 0) {
			rankedResult$isotope <- "M+0"
		} else if (iso == 1) {
			rankedResult$isotope <- "M+1"
		} else if (iso == 2) {
			rankedResult$isotope <- "M+2"
		} else if (iso == 3) {
			rankedResult$isotope <- "M+3"
		}
	# pseudoMSMS flag
		if (is.null(pseudoSpec)) {
		rankedResult$pseudoMSMS <- "FALSE"
		} else { rankedResult$pseudoMSMS <- "TRUE"
		}
	} else {
	output <- rankScore(result, specMatch)
	rankedResult <- output$rankedResult
	rankedSpec <- output$rankedSpecMatch
	}

	# Save output of ranked annotations -------
	# message("Saving results in '/Annotations' folder...")
	# set dataset name
	DatasetName <- paste(libs, "_", ESImode, "_", "RAMClustR", sep = "")
	write.csv(rankedResult[rankedResult$rank <= ncandidates,],
	          file = paste(mainDir, "/", subDir, "/", DatasetName, "_",
	                       round(fmz,3), "mz_", round(frt,3), "s_",
	                       "ranked_candidates.csv" , sep = ""), row.names = FALSE)

		# store result for high rank candidate
	if(exists('rankedResult') & !is.null(rankedResult)){
	global[i,"isotope"] <- as.character(rankedResult[1,"isotope"])
	global[i,"metabolite"] <- as.character(rankedResult[1,"metabolite"])
	global[i,"mz.metabolite"] <- rankedResult[1,"mz.metabolite"]
	global[i,"matched.mz"] <- rankedResult[1,"matched.mz"]
	global[i,"mz.error"] <- rankedResult[1,"mz.error"]
	global[i,"ion.type"] <- as.character(rankedResult[1,"ion.type"])
	global[i,"feature.type"] <- as.character(rankedResult[1,"feature.type"])
	global[i,"pseudoMSMS"] <- rankedResult[1,"pseudoMSMS"]
	global[i,"fraction"] <- as.character(rankedResult[1,"fraction"])
	global[i,"score"] <- rankedResult[1,"score"]
	} else {
	  global[i,c("metabolite", "feature.type", "ion.type", "isotope",
	             "mz.metabolite", "matched.mz", "mz.error", "pseudoMSMS",
	             "fraction", "score")] <- NA
	}
	# save image with matched spectra
	if (exists("rankedSpec") & !is.null(result)){
	  if (length(rankedSpec) == 1) {
	  plotCandidatesRC(fmz, frt, highCESpec, DatasetName, output, 1,
	                   DirPath = paste(mainDir,"/",subDir, "/", sep = ""))
	  }
	  if (length(rankedSpec) <= ncandidates) {
	      plots <- lapply(1:length(rankedSpec),
	                      function(x) plotCandidatesRC(fmz,
	                                                   frt,
	                                                   highCESpec,
	                                                   DatasetName,
	                                                   output, x,
	                                                   DirPath = paste(mainDir,
	                                                                   "/",
	                                                                   subDir,
	                                                                   "/",
	                                                                   sep = "")))
	  }
	  if (length(rankedSpec) > ncandidates) {
	    plots <- lapply(1:length(rankedSpec[1:ncandidates]),
	                    function(x) plotCandidatesRC(fmz,
	                                                 frt,
	                                                 highCESpec,
	                                                 DatasetName,
	                                                 output,
	                                                 x,
	                                                 DirPath = paste(mainDir,
	                                                                 "/",
	                                                                 subDir,
	                                                                 "/",
	                                                                 sep = "")))
	  }
	}	else NULL

	# save pseudo-MS/MS spectrum per feature
	if(exists("pseudoSpec") & !is.null(pseudoSpec) & length(pseudoSpec) > 0){
	  if (is.vector(pseudoSpec)) pseudoSpec <- as.data.frame(t(pseudoSpec))
	    pdf(file = paste(mainDir, "/", subDir, "/", "pseudoMSMS_",
	                     DatasetName, "_", round(fmz,3), "mz_", round(frt,3),
	                     "s", ".pdf", sep=""), width = 8, height = 5)

	plot(pseudoSpec[,1], pseudoSpec[,2], type = 'h',
	     xlim = c(50, max(pseudoSpec[,1]) + 100),
	     ylim = c(0, max(pseudoSpec[,2]) + max(pseudoSpec[,2])/1.5),
	     xlab = "m/z", ylab = "intensity (a.u.)", col = "black", lwd = 1,
	     main = paste("Pseudo-MS/MS Feature:", fmz, "m/z,", frt, "s"),
	     cex.main = 0.95, bty = "L", xaxs = "i", yaxs = "i")
	text(pseudoSpec[,1] - 10, pseudoSpec[,2],
	     as.character(round(pseudoSpec[,1],3)), pos = 4, cex = 0.8, srt = 45)
	dev.off()
	write.csv(
	  pseudoSpec, file = paste(mainDir, "/", subDir, "/", "pseudoMSMS_",
	                           DatasetName, "_", round(fmz,3), "mz_",
	                           round(frt,3), "s", ".csv", sep = ""),
	  row.names = FALSE)
	 } else NULL
 }
    # save global results table
   write.csv(
     global,
     file = paste(mainDir, "/", subDir, "/", "Global_Results", ".csv", sep = ""),
     row.names = FALSE)
	# save general options
	df <- data.frame(targetsTable_file = targetTable,
	                 RAMClusterObject = as.character(substitute(ramclustObj)),
	                 libraries = libs,
	                 ESImode = ESImode,
	                 RTfile = RTfile,
	                 checkIsotope = checkIsotope,
	                 matchWeight = matchWeight,
	                 tolerance = paste(tolerance, "ppm"),
	                 maxMZdiff = paste(maxMZdiff, "Da"),
	                 row.names = "Option")
	df <- as.data.frame(t(df))
	write.csv(df,
	          file = paste(mainDir, "/", subDir,
	                       "/", "General_options", ".csv", sep = ""))
message('Job done!')
}
