#' Obtain pseudo-MS/MS spectra for an LC-MS feature of interest.
#'
#' Function to obtain in-source MS and pseudo MS/MS spectra from a feature of
#' interest from All-ion fragmentation experiments (e.g. MSe, bbCID, AIF)
#'
#' @author  Goncalo Graca (Imperial College London)
#'
#' @param fmz The m/z for the feature of interest.
#' @param frt Retention time in seconds for the feature of interest.
#' @param xcmsF1 MSn object containing the LC-MS no-collision energy scans.
#' @param xcmsF2 MSn object containing the LC-MS all-ion fragmentation scans.
#' @param peaksF1 LC-MS picked peaks from xcmsF1 dataset using XCMS.
#' @param peaksF2 LC-MS picked peaks from xcmsF2 dataset using XCMS.
#' @param filetype The type of raw chromatogram imported: "mzML" or
#' "CDF" this is needed for scan frequency calculation.
#' @param nCE Number of collision energy MS2 functions contained in the data.
#' @param cthres1 Correlation threshold for the selection of in-source ions
#' related to the feature of interest.
#' @param cthres2 Correlation threshold for the selection of all-ion fragment
#' ions related to the feature of interest.
#' @param plotResults Logical argument indicate if the EICs and pseudo
#' in-source MS spectrum plots should be saved to disk.
#' @param savePseudoMSMS Logical argument indicate if the pseudo MS/MS should 
#' be saved to disk as .csv
#' @param scanfreq MS scanning frequency in seconds calculated as a the
#' difference between two consecutive low and high collision energy scans using
#' the information contained in xcmsF1 and xcmsF2 objects
#' @return A list containing several objects: insource, all MS1 peaks related to
#' the feature of interest; aif, all MS2 peaks related to the feature;
#' ms1_peaks, all MS1 peaks at the feature RT; ms2_peaks, all MS2 peaks at the
#' feature RT; ms2_eic, all EICs for the AIF features in the  RT window of the
#' feature of interest; mz_ms2, vector of m/z values for the MS2 ions in the
#' RT window of the feature of interest; feic, EIC of the feature of interest.
#' @export
getPseudoMSMS <- function(fmz, frt, xcmsF1, xcmsF2, peaksF1, peaksF2,
                          filetype = filetype, nCE = 1, cthres1 = 0.9,
                          cthres2 = 0.8, plotResults = TRUE, 
                          savePseudoMSMS = TRUE, SpName = "LCMS",
                          DirPath = paste(getwd(), "/", sep ="")){
  # create objects to store results from EIC correlations and peak-picking
  # improves object handling in other functions
  insource <- NULL
  aif <- NULL
  ms1_peaks <- NULL
  ms2_peaks <- NULL
  eic_aif <- NULL
  mz_ms2 <- NULL
  feic <- NULL

  # find feature
  delta <- 5
  fpeak <- chromPeaks(peaksF1, mz = fmz, ppm = delta)
  if(nrow(fpeak) == 0){
    while(nrow(fpeak) == 0){
      delta <- delta + 1
      fpeak <- chromPeaks(peaksF1, mz = fmz, ppm = delta)
    }
  }


  if (dim(fpeak)[1]==0) {
	message("feature not found in the peaks list")
  } else {
	if (dim(fpeak)[1]==1) i <- 1
	if (dim(fpeak)[1]>1) i <- which.min(abs(frt-fpeak[,"rt"]))

	# get rt, rtmin and rtmax
	rt <- fpeak[i,"rt"]
	rtmin <- fpeak[i,"rtmin"]
	rtmax <- fpeak[i,"rtmax"]
	mzmin <- fpeak[i,"mzmin"]
	mzmax <- fpeak[i,"mzmax"]

	# determine scan frequency in seconds------------------
	if(filetype == "mzML"){
	  if(nCE == 1){
	    spIdxF2 <- xcmsF2@featureData@data[1,"spIdx"]
	    spIdxF1 <- spIdxF2-1
	    idxF1 <- which(xcmsF1@featureData@data[,"spIdx"] == spIdxF1)
	    scanfreq <- xcmsF2@featureData@data[1,"retentionTime"] -
	      xcmsF1@featureData@data[idxF1,"retentionTime"]
	  }
	  if(nCE > 1) scanfreq <- xcmsF1@featureData@data[1,"retentionTime"] -
	      xcmsF2@featureData@data[1,"retentionTime"]
	 }
	if(filetype == "CDF"){
	  scanfreq <- xcmsF2@featureData@data[1,"retentionTime"] -
	    xcmsF1@featureData@data[1,"retentionTime"]
	}

	# function to correlate ms1 and/or ms2 EICs-----------------------
	eic_correlation <- function(a, b, scanfreq){
		fmaxRT <- rtime(a)[which.max(intensity(a))]
		tmaxRT <- rtime(b)[which.max(intensity(b))]
		inta <- intensity(a)
		intb <- intensity(b)
		if(!is.numeric(intb) | !is.numeric(inta)) {
		  c <- 0
		} else if(length(tmaxRT) == 0 | length(fmaxRT) == 0) {
		  c <- 0
		} else if(abs(tmaxRT - fmaxRT) > scanfreq) {
		  c <- 0
		} else if(length(na.omit(intb)) < 3) {
		  c <- 0
		} else if(length(inta) == length(intb)) {
		  c <- cor(inta, intb, use = "pairwise.complete.obs")
		} else if(length(inta) < length(intb)) {
		  len <- length(intb) - length(inta)
		  c <- cor(inta, intb[-(1:len)], use = "pairwise.complete.obs")
		} else if(length(inta) > length(intb)) {
		  len <- length(inta) - length(intb)
		   c <- cor(inta[-(1:len)], intb, use = "pairwise.complete.obs")
		}
    return(c)
	}

	# get in-source ions related with fmz----------------------------------------
	ms1_peaks <- xcms::chromPeaks(peaksF1)
	ms1_peaks <- ms1_peaks[which(ms1_peaks[,"rtmin"]<rt &
	                               ms1_peaks[,"rtmax"]>rt),]
	# to reduce the number of peaks and speed up the process...
	# ms1_peaks <- ms1_peaks[which(ms1_peaks[,"mz"] < fmz + 50),]
	# if (dim(ms1_peaks)[1] > 100)
	# ms1_peaks <- ms1_peaks[which(ms1_peaks[,"mz"] > fmz - 50),]
	npeaks <- dim(ms1_peaks)[1]

	# calculate EICs for all aparently coeluting peaks
	eic_ms1 <- lapply(1:npeaks,
	                  function(x) xcms::chromatogram(xcmsF1,
	                                           mz = c(ms1_peaks[x,"mz"] - 0.01,
	                                                  ms1_peaks[x,"mz"] + 0.01),
	                                           rt = c(rtmin - 1, rtmax + 1)))

	feic <- xcms::chromatogram(xcmsF1,
	                           mz = c(fmz - 0.01, fmz + 0.01),
	                           rt = c(rtmin - 1, rtmax + 1))

	c <- lapply(1:npeaks,
	            function(x) eic_correlation(feic[1,1],
	                                        eic_ms1[[x]][1,1], scanfreq))

	insource <- ms1_peaks[which(c > cthres1),]
	if (length(insource) > 12){
	  # sort features according to m/z:
		insource <- insource[order(insource[,"mz"], decreasing = FALSE),]
	} else NULL

	# get all MS2 peaks close to RT of MS1----------------------------------------
	# first need to get all MS2 peaks from peaksF2
	ms2_peaks <- xcms::chromPeaks(peaksF2)
	# then narrow down the selection to rtmin and rtmax of MS1 feature:
	ms2_peaks <- ms2_peaks[which(ms2_peaks[,"rtmin"] < rt &
	                               ms2_peaks[,"rtmax"] > rt),]
	mpeaks <- dim(ms2_peaks)[1]
	# if (mpeaks > 100) {
	# ms2_peaks <- ms2_peaks[which(ms2_peaks[,"mz"] < fmz + 50),]
	# } else NULL
	# to use at a later after matching peaks against the databases
	mz_ms2 <- ms2_peaks[,"mz"]
	mpeaks <- dim(ms2_peaks)[1]

	eic_aif <- lapply(1:mpeaks,
	                  function(x) xcms::chromatogram(
	                    xcmsF2, mz = c(ms2_peaks[x,"mz"]-0.01,
	                                   ms2_peaks[x,"mz"]+0.01),
	                    rt = c(rtmin-1,rtmax+1)))
	c2 <- lapply(1:mpeaks,
	             function(x) eic_correlation(feic[1,1],
	                                         eic_aif[[x]][1,1], scanfreq))

	aif <- ms2_peaks[which(c2 > cthres2),]

	if(length(aif) == 0) aif <- NULL

	# plot EICs and pseudo-MS using ggplot2 and gridExtra-------------------------

	if(plotResults & !is.null(aif)) {
		if (length(insource) > 12){
			is_ions <- match(insource[,"mz"], ms1_peaks[,"mz"])
		} else is_ions <- match(insource["mz"], ms1_peaks[,"mz"])

    df1 <- lapply(is_ions,
                  function(x) data.frame(
                    intensity = xcms::intensity(eic_ms1[[x]][1,1]),
                    rt = xcms::rtime(eic_ms1[[x]][1,1]), mz = as.character(
                      rep(paste(round(ms1_peaks[x,"mz"],3),"m/z")))))
    df1 <- do.call("rbind", df1)

    p1 <- ggplot2::ggplot(df1[!is.na(df1$intensity),],
                          ggplot2::aes(x = rt,
                                       y = intensity, colour = as.factor(mz))) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::labs(x = "RT (s)", y = "Intensity (a.u.)",
                    colour = paste("Correlation >", cthres1)) +
      ggplot2::ggtitle(paste("Correlated EICs:", round(fmz,3),
                             "m/z",round(frt),"s")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    if (length(insource) > 12){
      df2 <- as.data.frame(insource[,c("mz", "into")])
    } else df2 <- data.frame(mz = insource["mz"], into = insource["into"])

    p2 <- ggplot2::ggplot(df2,
                          ggplot2::aes(x = mz, y = into,
                                       label = round(mz, 3))) +
      ggplot2::geom_segment(ggplot2::aes(xend = mz, yend=0),
                            color="red", lwd=0.5) +
      ggplot2::geom_text(size=3, angle=45, hjust=0, vjust=0) +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle(paste("Pseudo-MS:", round(fmz,3),"m/z",round(frt),"s")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::ylim(0, max(df2[,2]) + 0.1*max(df2[,2])) +
      ggplot2::xlim(min(df2[,1])-50, max(df2[,1])+50) +
      ggplot2::labs(x = "m/z", y = "Intensity (a.u.)")

    if (length(aif) > 12){
      aif_ions <- match(aif[,"mz"], ms2_peaks[,"mz"])
    } else aif_ions <- match(aif["mz"], ms2_peaks[,"mz"])

    df3 <- lapply(aif_ions,
                  function(x) data.frame(
                    intensity = xcms::intensity(eic_aif[[x]][1,1]),
                    rt = xcms::rtime(eic_aif[[x]][1,1]), mz = as.character(
                      rep(paste(round(ms2_peaks[x,"mz"],3),"m/z")))))
    df3 <- do.call("rbind", df3)

    p3 <- ggplot2::ggplot(df3[!is.na(df3$intensity),],
                          ggplot2::aes(x = rt, y = intensity, colour = mz)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::labs(x = "RT (s)", y = "Intensity (a.u.)",
           colour = paste("Correlation >", cthres2)) +
      ggplot2::ggtitle(paste("Fragments EICs:", round(fmz,3),
                             "m/z",round(frt),"s")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    if (length(aif) > 12){
      df4 <- as.data.frame(aif[,c("mz", "into")])
    } else df4 <- data.frame(mz = aif["mz"], into = aif["into"])

    p4 <- ggplot2::ggplot(df4,
                          ggplot2::aes(x=mz, y=into, label = round(mz, 3))) +
      ggplot2::geom_segment(ggplot2::aes(xend = mz, yend=0),
                            color="red", lwd=0.5) +
      ggplot2::geom_text(size=3, angle=45, hjust=0, vjust=0) +
      ggplot2::ggtitle(paste("Pseudo-MS/MS:", round(fmz,3),
                             "m/z",round(frt),"s")) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::ylim(0, max(df4[,2]) + 0.1*max(df4[,2])) +
      ggplot2::xlim(min(df4[,1])-50, max(df4[,1])+50) +
      ggplot2::labs(x = "m/z", y = "Intensity (a.u.)")

    pdf(paste(DirPath, SpName, "_", "pseudoMS_AIF_", round(fmz,3),
              "mz_", round(frt), "s", ".pdf", sep = ""), height = 8, width = 12)
    gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)
    dev.off()
    
    if(savePseudoMSMS){
      if(!is.null(aif)) {
        write.csv(df4, paste(DirPath, SpName, "_", "pseudoMS_AIF_", 
                             round(fmz, 3), "mz_",round(frt),"s", 
                             ".csv", sep = ""))
      }
    }

	} else {

    if (length(insource) > 12){
		is_ions <- match(insource[,"mz"], ms1_peaks[,"mz"])
    } else is_ions <- match(insource["mz"], ms1_peaks[,"mz"])

    df1 <- lapply(is_ions,
                  function(x) data.frame(
                    intensity = xcms::intensity(eic_ms1[[x]][1,1]),
                    rt = xcms::rtime(eic_ms1[[x]][1,1]),
                    mz = as.character(
                      rep(paste(round(ms1_peaks[x,"mz"],3),"m/z")))))
    df1 <- do.call("rbind", df1)

    p1 <- ggplot2::ggplot(df1[!is.na(df1$intensity),],
                          ggplot2::aes(x = rt, y = intensity,
                                       colour = as.factor(mz))) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::labs(x = "RT (s)", y = "Intensity (a.u.)",
           colour = paste("Correlation >", cthres1)) +
      ggplot2::ggtitle(paste("Correlated EICs:",
                             round(fmz,3),"m/z",round(frt),"s")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

	if (length(insource) > 12){
      df2 <- as.data.frame(insource[,c("mz", "into")])
    } else df2 <- data.frame(mz = insource["mz"], into = insource["into"])

    p2 <- ggplot2::ggplot(df2,
                          ggplot2::aes(x = mz, y = into, label = round(mz,3))) +
      ggplot2::geom_segment(ggplot2::aes(xend = mz, yend=0),
                            color="red", lwd=0.5) +
      ggplot2::geom_text(size=3, angle=45, hjust=0, vjust=0) +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle(paste("Pseudo-MS:", round(fmz,3),"m/z",round(frt),"s")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::ylim(0, max(df2[,2]) + 0.1*max(df2[,2])) +
      ggplot2::xlim(min(df2[,1])-50, max(df2[,1])+50) +
      ggplot2::labs(x = "m/z", y = "Intensity (a.u.)")

    pdf(paste(DirPath, SpName, "_", "in-source_", round(fmz, 3),
              "mz_",round(frt),"s",".pdf", sep = ""), height = 8, width = 12)
    gridExtra::grid.arrange(p1, p2, nrow = 1)
    dev.off()
	}
  }
  result <- list(insource = insource,
                 aif = aif,
                 ms1 = ms1_peaks,
                 ms2 = ms2_peaks,
                 ms2_eic = eic_aif,
                 mz_ms2 = mz_ms2,
                 feic = feic)
  return(result)
}
