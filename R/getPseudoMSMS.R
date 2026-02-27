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
#' Should be set to NULL to obtain only the in-source fragmentation (ISF) 
#' pseudo-MS/MS.
#' @param peaksF1 LC-MS picked peaks from xcmsF1 dataset using XCMS.
#' @param peaksF2 LC-MS picked peaks from xcmsF2 dataset using XCMS. 
#' Should be set to NULL to obtain only the in-source fragmentation (ISF) 
#' pseudo-MS/MS.
#' @param cthres1 Correlation threshold for the selection of in-source ions
#' related to the feature of interest.
#' @param cthres2 Correlation threshold for the selection of all-ion fragment
#' ions related to the feature of interest.
#' @return A list containing several objects: insource, all MS1 peaks related 
#' to the feature of interest; aif, all MS2 peaks related to the feature;
#' ms1_peaks, all MS1 peaks at the feature RT; ms2_peaks, all MS2 peaks at the
#' feature RT; ms2_eic, all EICs for the AIF features in the  RT window of the
#' feature of interest; mz_ms2, vector of m/z values for the MS2 ions in the
#' RT window of the feature of interest; feic, EIC of the feature of interest;
#' feic_aif, the EICs of all MS2 ions correlated with the feature of interest.
#' If xcmsF2 is set to NULL the in-source pseudo-MS/MS spectrum will be saved 
#' instead of the AIF pseudo-MS/MS and similarly the EICs from the MS1 ions
#' correlated with the feature of interest.
#' @examples
#' # obtain the pseudo-MS/MS of one feature from the 
#' # MS1 scans (in-source fragments)
#' # read the example LC-MS data from the msdata package
#' library(msdata)
#' filePath <- system.file("TripleTOF-SWATH", "PestMix1_SWATH.mzML", 
#' package="msdata")
#' xcmsF1 <- MSnbase::readMSData(filePath, msLevel.=1, mode="onDisk")
#' # perform peak-picking
#' cwp <- CentWaveParam(snthresh=5, noise=100, ppm=10, peakwidth=c(3, 30))
#' peaksF1 <- xcms::findChromPeaks(xcmsF1, msLevel = 1L, param = cwp)
#' # feature m/z and RT
#' fmz <- 304.1124
#' frt <- 423.945
#' # obtain the pseudo-MS/MS from MS1 scans (in-source fragments spectrum):
#' getPseudoMSMS(fmz, frt, xcmsF1, xcmsF2=NULL, peaksF1, peaksF2=NULL, 
#' cthres1=0.9, cthres2=0.8)
#' @importFrom xcms chromPeaks chromatogram 
#' @importFrom MSnbase compareChromatograms
#' @importFrom ProtGenerics intensity rtime polarity
#' @export
getPseudoMSMS <- function(fmz, frt, xcmsF1, xcmsF2, peaksF1, peaksF2, 
						  cthres1=0.9, cthres2=0.8){
    ## create objects to store results from EIC correlations and peak-picking--
    # improves object handling in other functions
    insource <- NULL
    aif <- NULL
    ms1_peaks <- NULL
    ms2_peaks <- NULL
    eic_ms1 <- NULL
    eic_aif <- NULL
    mz_ms2 <- NULL
    feic <- NULL
    targetRT <- NULL
    ## find feature peak in peaksF1 object-------------------------------------
    # increase tolerance by 1 from 5 ppm until the feature is found
    delta <- 5
    fpeak <- chromPeaks(peaksF1, mz=fmz, ppm=delta)
    if(nrow(fpeak) == 0){
        while(nrow(fpeak) == 0 & delta < 30){
            delta <- delta + 1
            fpeak <- chromPeaks(peaksF1, mz=fmz, ppm=delta)
        }
    }
    
    if(nrow(fpeak) == 0) message("feature not found in the peaks list")
    
    if(is.matrix(fpeak) & nrow(fpeak) == 1) targetRT <- fpeak[,"rt"]
    if(nrow(fpeak) > 1) {
        idx <- which.min(abs(frt-fpeak[,"rt"]))
        # get rt, rtmin and rtmax
        fpeak <- fpeak[idx,]
        fpeak <- as.matrix(t(fpeak))
        targetRT <- fpeak[,"rt"]
    }

    ## get in-source ions related with fmz-------------------------------------
    if(!is.null(targetRT)){
        isf <- getISF(xcmsF1, peaksF1, fpeak, targetRT, cthres1)
        ms1_peaks <- isf$ms1_peaks 
        feic <- isf$feic 
        eic_ms1 <- isf$eic_ms1
        insource <- isf$insource
    }
    ## get all MS2 peaks close to RT of MS1 and pseudo MS/MS spectrum----------
    if(!is.null(xcmsF2) & !is.null(peaksF2)){
        AIF <- getAIF(fmz, frt, xcmsF2, peaksF2, targetRT, feic, cthres2)
        aif <- AIF$aif
        ms2_peaks <- AIF$ms2_peaks
        eic_aif <- AIF$eic_aif
        mz_ms2 <- AIF$mz_ms2
    } else eic_aif <- eic_ms1 # stores ISF EICs instead for later plotting

    ## return pseudo-MS/MS object
    result <- list(insource=insource, aif=aif, ms1=ms1_peaks, ms2=ms2_peaks,
                    ms2_eic=eic_aif, mz_ms2=mz_ms2, feic=feic)
    return(result)
}


## Helper functions------------------------------------------------------------

## Get in-source fragments and other feature-related objects-------------------
getISF <- function(xcmsF1, peaksF1, fpeak, targetRT, cthres1){
    ms1_peaks <- xcms::chromPeaks(peaksF1)
    ms1_peaks <- ms1_peaks[which(ms1_peaks[,"rtmin"] < targetRT &
                                    ms1_peaks[,"rtmax"] > targetRT),]
    npeaks <- nrow(ms1_peaks)

    # get EIC of the feature of interest
    feic <- xcms::chromatogram(xcmsF1,
                                mz=fpeak[,c("mzmin","mzmax")],
                                rt=fpeak[,c("rtmin","rtmax")])
    
    # calculate EICs for all apparently co-eluting peaks
    eic_ms1 <- xcms::chromatogram(xcmsF1, mz = ms1_peaks[,c("mzmin", "mzmax")],
                                    rt = ms1_peaks[,c("rtmin","rtmax")])
    
    # correlate feature EIC with co-eluting peak EICs
    c <- MSnbase::compareChromatograms(eic_ms1,
                                        feic,
                                        ALIGNFUNARGS=list(method="closest"))
    
    insource <- ms1_peaks[which(c > cthres1),]
    if(length(insource) > 12){
        # sort features according to m/z:
        insource <- insource[order(insource[,"mz"], decreasing=FALSE),]
    } else NULL
    isf <- list(ms1_peaks=ms1_peaks, feic=feic, 
                eic_ms1=eic_ms1, insource=insource)
    return(isf)
}

## Function to obtain relevant AIF fragments and other AIF-related objects-----
getAIF <- function(fmz, frt, xcmsF2, peaksF2, targetRT, feic, cthres2){
    # get all MS2 peaks close to RT of MS1
    # first need to get all MS2 peaks from peaksF2
    ms2_peaks <- xcms::chromPeaks(peaksF2)
    # then narrow down the selection to rtmin and rtmax of MS1 feature:
    ms2_peaks <- ms2_peaks[which(ms2_peaks[,"rtmin"] < targetRT &
                                ms2_peaks[,"rtmax"] > targetRT),]
    mz_ms2 <- ms2_peaks[,"mz"]
    mpeaks <- nrow(ms2_peaks)
    
    eic_aif <- xcms::chromatogram(xcmsF2, 
                                    mz = ms2_peaks[,c("mzmin","mzmax")],
                                    rt = ms2_peaks[,c("rtmin","rtmax")],
                                    msLevel = 2)
    c2 <- MSnbase::compareChromatograms(eic_aif,
                                        feic,
                                        ALIGNFUNARGS=list(method="closest"))
    aif <- ms2_peaks[which(c2 > cthres2),]
    if(length(aif) == 0) aif <- NULL
    
    # remove ions which are much greater than fmz by 5 amu
    if(!is.null(aif)){
        aif <- aif[which(aif[,"mz"] < fmz + 5),]
    }
    result <- list(aif=aif, 
                    ms2_peaks=ms2_peaks,
                    eic_aif=eic_aif,
                    mz_ms2=mz_ms2)
    return(result)
}