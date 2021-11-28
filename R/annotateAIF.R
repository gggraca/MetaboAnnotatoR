#' Function to annotate features using LC-MS AIF raw chromatograms.
#'
#' Search possible matches of a feature in Lipid fragments and other small
#' molecule libraries, using 'raw' LC-MS single chromatogram information
#' already loaded and peak-picked.
#'
#' @author Goncalo Graca (Imperial College London)
#'
#' @param targetTable A .csv file containing the list of features to annotate
#' and the name of the files containing the raw data.
#' @param ESImode Ionization mode: 'POS' for positive (default)
#' or 'NEG' for negative ionisation modes
#' @param libs Fragment libraries to use: 'Lipids' (default) or 'Metabolites'
#' for other small molecules.
#' @param RTfile Optional csv file with Lipid/metabolites classes Retention
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
#' @param ncandidates Maximum number of candidates to plot and store.
#' @return For each feature in the targeTable the will return a ranked list of
#' annotations, a plot of EICs and pseudo-MS/MS spectrum for the matched ions,
#' a plot of pseudo-MS/MS and pseudo-MS spectra for the each feature,
#' a targeTable annotated with rank 1 annotations and a table with the options
#' used for the function.
#' @examples
#' # get the some example human serum LC-MS data and feaure list to annotate:
#' getDemoData()
#' # run the annotation using the lipid libraries:
#' annotateAIF(targetTable = "targetTable.csv", filetype = "mzML", libs = "Lipids",
#' ESImode = "POS", RTfile = "none", nCE = 1, corThresh = 0.7, checkIsotope = TRUE)
#' @export
annotateAIF <- function(targetTable = NULL,
                        filetype = "mzML",
                        libs = "Lipids",
                        ESImode = "POS",
                        RTfile = "none",
                        nCE = 1,
                        corThresh = 0.8,
					              checkIsotope = TRUE,
					              tolerance = 25,
					              maxMZdiff = 0.01,
					              matchWeight = 0.5,
					              ncandidates = 5){
# Read XCMS peak-picking options
if (file.exists("XCMS_options.csv")) {
  xcmsOptions <- read.csv("XCMS_options.csv")
} else {
  message("No XCMS options file found in the working directory")
  answer <- readline(prompt = "Use default options (Yes/No)? [y/n]")
  if(answer == "y") {
    xcmsOptionsPath <- system.file("XCMS_options.csv",
                                   package = "MetaboAnnotatoR")
    xcmsOptions <- read.csv(xcmsOptionsPath)
  }
  if(answer == "n"){
    xcmsOptionsPath <- system.file("XCMS_options.csv",
                                   package = "MetaboAnnotatoR")
    file.copy(from = xcmsOptionsPath, to = getwd())
    message("Default  XCMS options file saved in the working directory")
    stop("Please edit the file and re-run.")
  }
}

# Read targets table
  if(is.null(targetTable)){
    message("Targets table not found")
    targetTablePath <- system.file("targetTable.csv",
                                   package = "MetaboAnnotatoR")
    file.copy(from = targetTablePath, to = getwd())
    message("Default targetTable file saved in the working directory.")
    message("Please edit the file and replace it in the function argument.")
    stop()
    } else {
  targets <- read.csv(targetTable)
  }
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

# evaluate if one or more samples are to be read-----------
if(length(unique(targets[,3])) == 1){
	message("Reading data...")
  if(filetype == "mzML"){
	  data.path <- paste("./", targets[1,3], ".mzML", sep = "")
	  # separate the two MS functions
	  xcmsF1 <- MSnbase::readMSData(data.path, msLevel. = 1, mode = "onDisk")
	  xcmsF2 <- MSnbase::readMSData(data.path, msLevel. = 2, mode = "onDisk")
	  if(nCE > 1){
	    maxCE <- max(xcmsF2@featureData@data$collisionEnergy)
	    highCEscans <- which(xcmsF2@featureData@data$collisionEnergy == maxCE)
	    xcmsF2@featureData@data <- xcmsF2@featureData@data[highCEscans,]
	  }
	  # this is necessary to enable XCMS to read from MS2 scans
	  xcmsF2@featureData@data$msLevel <- 1
	}
  if(filetype == "CDF"){
	  data.path <- paste("./", targets[1,3], "01.CDF", sep = "")
	  data.path2 <- paste("./", targets[1,3], "02.CDF", sep = "")
	  # read the two MS functions
	  xcmsF1 <- MSnbase::readMSData(data.path, mode = "onDisk")
	  xcmsF2 <- MSnbase::readMSData(data.path2, mode = "onDisk")
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
		  xcmsF1 <- MSnbase::readMSData(data.path, msLevel. = 1, mode = "onDisk")
		  xcmsF2 <- MSnbase::readMSData(data.path, msLevel. = 2, mode = "onDisk")
		  xcmsF2@featureData@data$msLevel <- 1
	}
	if(filetype == "CDF"){
		data.path <- paste("./", targets[i,3], "01.CDF", sep = "")
		data.path2 <- paste("./", targets[i,3], "02.CDF", sep = "")
	# read the two MS functions
		xcmsF1 <- MSnbase::readMSData(data.path, mode = "onDisk")
		xcmsF2 <- MSnbase::readMSData(data.path2, mode = "onDisk")
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
	                         filetype = filetype, nCE = 1, cthres1 = corThresh,
	                         cthres2 = corThresh, plotResults = TRUE, 
	                         savePseudoMSMS = TRUE, SpName, DirPath)
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
	message("Mathing candidate(s) fragments to pseudo-MS/MS and highCE spectra...")
	output <- mapply(compFrag,
	                 candidates,
	                 lapply(as.numeric(names(candidates)), function(x) lib[[x]]),
	                 MoreArgs = list(fmz, frt, iso, highCESpec,
	                                 pseudoSpec,
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
message("Job done!")
}
