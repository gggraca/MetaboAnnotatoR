#' Function annotateAIF for raw LC-MS chromatogram
#' Search possible matches of a feature in Lipid fragments and other small
#' molecule libraries
#' using 'raw' LC-MS single chromatogram information already loaded and
#' peak-picked returns a ranked list of candidates for the feature annotation
#' as .csv file
#' Arguments:
#' targetTable: a .csv file containing the list of features to annotate and the
#' name of the files containing the raw data
#' ESImode: ionization mode: 'POS' for positive (default) or 'NEG' for negative
#' libs: fragments library: 'Lipids' (default) or 'Metabolites' for other
#' small molecules
#' RTfile: .csv file with Lipid/metabolites classes RT information (optional)
#' corThresh: Pearson correlation coefficient for EIC correlation
#' checkIsotope: whether or not to check the isotope type;
#' default is to check (TRUE)
#' tolerance: tolerance in ppm for the candidate search
#' maxMZdiff: maximum m/z difference between candidate fragments and
#' pseudo-MS/MS or AIF ions
#' Goncalo Graca & Yuheng Cai (Imperial College London)
#' @export
annotateAIF <- function(targetTable,
                        filetype = "mzML",
                        libs = "Lipids",
                        ESImode = "POS",
                        RTfile = "none",
                        corThresh = 0.8,
					              checkIsotope = TRUE,
					              tolerance = 25,
					              maxMZdiff = 0.01,
					              matchWeight = 0.5,
					              ncandidates = 5){
# Read XCMS peak-picking options
xcmsOptions <- read.csv("XCMS_options.csv")

# Read targets table
targets <- read.csv(targetTable, header = TRUE)

# create table to store global results
global <- targets
global[,c("metabolite", "feature.type", "ion.type", "isotope", "mz.metabolite",
"matched.mz", "mz.error", "pseudoMSMS", "fraction", "score")] <- NA
# Get RTs from file --------------------
	if(RTfile == "none") {
		message("No RT information provided...")
		RTs <- "none"
		} else {
		message("Reading RT information...")
		RTs <- read.csv(RTfile, header = TRUE)
		}

# Create directory to store the results-------------
mainDir <- "./Annotations"
Date <- Sys.Date()
Time <- format(Sys.time(), "%X")
Time <- gsub(":", "_", Time)
subDir <- paste(libs, "_", ESImode, "_", "AIF", "_", Date,"_", Time, sep="")
dir.create(file.path(mainDir), showWarnings = FALSE)
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
DirPath <- paste(mainDir, "/", subDir, "/", sep = "")

# load libraries --------------------
message("Importing libraries...")
libfiles <- list.files(path = paste("./Libraries/",libs,"/", ESImode, sep = ""),
                       full.names = TRUE)
# check.names=FALSE to use the original header names in the annotations:
lib <- sapply(libfiles, read.csv, header = TRUE, sep=",", check.names = FALSE)

# evaluate if one or more samples are to be read-----------
if(length(unique(targets[,3])) == 1){
	message("Reading data...")
  if(filetype == "mzML"){
	  data.path <- paste("./", targets[1,3], ".mzML", sep = "")
	  # separate the two MS functions
	  xcmsF1 <- xcms::readMSData(data.path, msLevel. = 1, mode = "onDisk")
	  xcmsF2 <- xcms::readMSData(data.path, msLevel. = 2, mode = "onDisk")
	  # this is necessary to enable xcms to read from MS2 scans
	  xcmsF2@featureData@data$msLevel <- 1
	}
  if(filetype == "CDF"){
	  data.path <- paste("./", targets[1,3], "01.CDF", sep = "")
	  data.path2 <- paste("./", targets[1,3], "02.CDF", sep = "")
	  # read the two MS functions
	  xcmsF1 <- xcms::readMSData(data.path, mode = "onDisk")
	  xcmsF2 <- xcms::readMSData(data.path2, mode = "onDisk")
	}
  # try to use same centwave parameters for both no- and high-collision scans
  cwp <- xcms::CentWaveParam(snthresh = xcmsOptions[3,2],
                        noise = xcmsOptions[2,2],
                        ppm = xcmsOptions[1,2],
                        peakwidth = as.numeric(xcmsOptions[4,2:3]),
                        prefilter = as.numeric(xcmsOptions[5,2:3]))
  peaksF1 <- xcms::findChromPeaks(xcmsF1, param = cwp)
  peaksF2 <- xcms::findChromPeaks(xcmsF2, param = cwp)
} else NULL


for (i in 1:dim(targets)[1]){
	message(paste("####### Processing feature", i, "of",
	              dim(targets)[1]), " ########")
	if(length(unique(targets[,3]))>1){
		message("Reading data...")
	  if(filetype == "mzML"){
		  data.path <- paste("./", targets[i,3], ".mzML", sep = "")
	  # separate the two MS functions
		  xcmsF1 <- xcms::readMSData(data.path, msLevel. = 1, mode = "onDisk")
		  xcmsF2 <- xcms::readMSData(data.path, msLevel. = 2, mode = "onDisk")
		  xcmsF2@featureData@data$msLevel <- 1
	}
	if(filetype == "CDF"){
		data.path <- paste("./", targets[i,3], "01.CDF", sep = "")
		data.path2 <- paste("./", targets[i,3], "02.CDF", sep = "")
	# read the two MS functions
		xcmsF1 <- xcms::readMSData(data.path, mode = "onDisk")
		xcmsF2 <- xcms::readMSData(data.path2, mode = "onDisk")
	}
	  # try to use same centwave parameters for both no- and high-collision scans
	  cwp <- xcms::CentWaveParam(snthresh = xcmsOptions[3,2],
	                       noise = xcmsOptions[2,2],
	                       ppm = xcmsOptions[1,2],
	                       peakwidth=as.numeric(xcmsOptions[4,2:3]),
	                       prefilter=as.numeric(xcmsOptions[5,2:3]))
	  peaksF1 <- xcms::findChromPeaks(xcmsF1, param = cwp)
	  peaksF2 <- xcms::findChromPeaks(xcmsF2, param = cwp)
	} else NULL

	fmz <- targets[i,1]
	frt <- targets[i,2]

	# get sample name to save files into ---------------
	Sp <- strsplit(data.path, "/")
	Sp.idx <- length(Sp[[1]])
	SpName <- Sp[[1]][Sp.idx]
	if(filetype == "mzML") SpName <- gsub(".mzML", "", SpName)
	if(filetype == "CDF") SpName <- gsub("01.CDF", "", SpName)

	# get MS spectra at feature RT --------------------
	message("Obtaining pseudo-MS/MS spectrum...")
	try(
	  specs <- getPseudoMSMS(fmz, frt, xcmsF1, xcmsF2, peaksF1, peaksF2,
	                         filetype = filetype, CE = 1, cthres1 = corThresh,
	                         cthres2 = corThresh, plotResults = TRUE, SpName,
	                         DirPath)
	)

	if(exists("specs")){
	  highCESpec <- specs$ms2
	  pseudoSpec <- specs$aif
	  inSourceSpec <- specs$insource
	  ms2eic <- specs$ms2_eic
	  feic <- specs$feic
	} else next

	# Isotope check --------------------
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

	# Search Libraries --------------------
	if(is.null(pseudoSpec) & is.null(highCESpec)) { next
	} else {
	message("Searching candidates...")
	candidates <- searchLib(lib, libfiles, fmz - iso * 1.0034,
	                        frt, tolerance = tolerance, RTs, inSourceSpec)
	}

	# Compare fragments between Library candidates and
	# high-collision-energy / pseudo-MS/MS spectra --------------------
	if(is.null(pseudoSpec) & is.null(highCESpec) |
	   length(unlist(candidates)) == 0) {
	  result <- NULL
	} else {
	message("Comparing pseudo MS/MS and high collision energy MS
	        with candidate(s) fragments...")
	output <- mapply(compFrag,
	                 candidates,
	                 lapply(as.numeric(names(candidates)), function(x) lib[[x]]),
	                 MoreArgs = list(fmz, frt, iso, highCESpec,
	                                 pseudoSpec, tolerance = tolerance,
	                                 maxMZdiff = maxMZdiff,
	                                 matchWeight = matchWeight),
	                 SIMPLIFY = FALSE)
	result <- do.call(rbind,lapply(output, "[[", 1))
	specMatch <- unlist(lapply(output, "[[",2), recursive = FALSE)
	specMatch <- specMatch[!(specMatch) == "NULL"]
	}
	# Score ranking --------------------
	if(is.null(result)) {
	rankedResult <- targets[i, 1:2]
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
	output <- rankScore(result,specMatch)
	rankedResult <- output$rankedResult
	rankedSpec <- output$rankedSpecMatch
	}
	# Save output of ranked annotations: ranks 1:5 -------
	# message("Saving results in '/Annotations' folder...")
	if(is.null(result)) {
	write.csv(rankedResult,
	          file = paste(mainDir, "/",
	                       subDir,"/",
	                       SpName,"_", round(fmz,3),"mz_",
	                       round(frt,3),"s_",
	                       "ranked_candidates.csv",sep = ""),
	          row.names = FALSE)
	} else {
	write.csv(rankedResult[rankedResult$rank <= ncandidates,],
	          file = paste(mainDir, "/", subDir,"/", SpName,"_",
	                       round(fmz,3),"mz_", round(frt,3),"s_",
	                       "ranked_candidates.csv",sep = ""), row.names = FALSE)
	}
	# store result for high rank candidate
	if(exists("rankedResult") & !is.null(rankedResult)){
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
	# save image with matched spectra and EICs
	# must update the folder to save images in...
	# can be updated to work with single feature
	if (exists('rankedSpec') & !is.null(result)){
	  if (length(rankedSpec) == 1) {
	    plotCandidatesAIF(fmz, frt, highCESpec,
	                      ms2eic, SpName, output, 1, DirPath)
	  }
	  if (length(rankedSpec) <= ncandidates) {
	    #try(
	      plots <- lapply(1:length(rankedSpec),
	                      function(x) plotCandidatesAIF(fmz, frt,
	                                                    highCESpec,
	                                                    ms2eic, SpName,
	                                                    output, x, DirPath))
	    #)
	  }
	  if (length(rankedSpec) > ncandidates) {
	    #try(
	    plots <- lapply(1:length(rankedSpec[1:ncandidates]),
	                    function(x) plotCandidatesAIF(fmz, frt,
	                                                  highCESpec,
	                                                  ms2eic, SpName,
	                                                  output, x, DirPath))
	    #)
	  }
	}	else NULL
  }
    # save global results table
   write.csv(global,
             file = paste(mainDir, "/", subDir, "/", "Global_Results",
                          ".csv", sep = ""), row.names = FALSE)
	# save general options
	df <- data.frame(targetsTable_file = targetTable, libraries = libs,
	                 ESImode = ESImode, RTfile = RTfile, corThresh = corThresh,
	                 checkIsotope = checkIsotope, matchWeight = matchWeight,
	                 tolerance = paste(tolerance, "ppm"),
	                 maxMZdiff = paste(maxMZdiff, "Da"), row.names = "Option")
	df <- as.data.frame(t(df))
	write.csv(df, file = paste(mainDir, "/", subDir, "/",
	                           "General_options", ".csv", sep = ""))
	# save xcms options
	write.csv(xcmsOptions,
	          file = paste(mainDir, "/", subDir, "/", "XCMS_options",
	                       ".csv", sep = ""), row.names = FALSE)
message('Job done!')
}
