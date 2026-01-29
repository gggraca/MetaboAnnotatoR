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
#' @param filetype The type of raw chromatogram imported: "mzML" or
#' "CDF" this is needed for scan frequency calculation.
#' @param nCE Number of collision energy MS2 functions contained in the data.
#' @param cthres1 Correlation threshold for the selection of in-source ions
#' related to the feature of interest.
#' @param cthres2 Correlation threshold for the selection of all-ion fragment
#' ions related to the feature of interest.
#' @param savePlotResults Logical argument indicate if the EICs and pseudo
#' in-source MS spectrum plots should be saved to disk.
#' @param savePseudoMSMS Logical argument indicate if the pseudo MS/MS should 
#' be saved to disk as .csv.
#' @param ExpName Name of the experimental technique, commonly LC-MS (default). 
#' Other names are possible, e.g. CE-MS, LC-AIF-MS, etc..
#' @param DirPath Path to the folder where plots will be saved. 
#' The default folder is the working directory.
#' @return A list containing several objects: insource, all MS1 peaks related 
#' to the feature of interest; aif, all MS2 peaks related to the feature;
#' ms1_peaks, all MS1 peaks at the feature RT; ms2_peaks, all MS2 peaks at the
#' feature RT; ms2_eic, all EICs for the AIF features in the  RT window of the
#' feature of interest; mz_ms2, vector of m/z values for the MS2 ions in the
#' RT window of the feature of interest; feic, EIC of the feature of interest.
#' If xcmsF2 are both set to NULL and savePseudoMSMS is TRUE the in-source 
#' pseudo-MS/MS spectrum will be saved instead.
#' @examples
#' # obtain the pseudo-MS/MS of one feature from the 
#' MS1 scans (in-source fragments)
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
#' filetype="mzML", nCE=1, cthres1=0.9, cthres2=0.8, savePlotResults=FALSE, 
#' savePseudoMSMS=FALSE, ExpName="LCMS", DirPath=paste(getwd(), "/", sep =""))
#' @importFrom xcms chromPeaks chromatogram 
#' @importFrom MSnbase compareChromatograms
#' @importFrom ProtGenerics intensity rtime polarity
#' @export
getPseudoMSMS <- function(fmz, frt, xcmsF1, xcmsF2, peaksF1, peaksF2,
                            filetype=filetype, nCE=1, cthres1=0.9,
                            cthres2=0.8, savePlotResults=TRUE, 
                            savePseudoMSMS=TRUE, ExpName="LCMS",
                            DirPath=paste(getwd(), "/", sep ="")){
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

    ## plot in-source fragmentation (ISF) and AIF pseudo-MS/MS spectra---------
    if(savePlotResults & !is.null(aif)) {
        plotPseudoMSMS(fmz, frt,
                        insource, ms1_peaks, eic_ms1, 
                        aif, ms2_peaks, eic_aif, cthres1, cthres2,
                        DirPath, ExpName)
    }

    ## save AIF pseudo-MS/MS spectrum as .csv and .mgf-------------------------
    if(savePseudoMSMS & !is.null(aif)){
        savePseudoMSMS(fmz, frt, xcmsF1, pseudo=aif, specType="AIF", DirPath, 
                        ExpName)
    } else if(savePseudoMSMS & is.null(aif) & !is.null(insource)){
        plotISF(fmz, frt, insource, ms1_peaks, eic_ms1, cthres1, DirPath, 
                ExpName)
    }
    
    ## save ISF pseudo-MS/MS spectrum as .csv and .mgf-------------------------
    if(savePseudoMSMS & is.null(aif) & !is.null(insource)){
        savePseudoMSMS(fmz, frt, xcmsF1, pseudo=insource, specType="ISF", 
                        DirPath, ExpName)
    }

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

## Plot AIF and ISF pseudo-MS/MS results to pdf--------------------------------
plotPseudoMSMS <- function(fmz, frt, 
                            insource, ms1_peaks, eic_ms1, 
                            aif, ms2_peaks, eic_aif, cthres1, cthres2,
                            DirPath, ExpName){
    if(length(insource) > 12){
        is_ions <- match(insource[,"mz"], ms1_peaks[,"mz"])
    } else is_ions <- match(insource["mz"], ms1_peaks[,"mz"])
    
    df1 <- lapply(is_ions,
                    function(x) data.frame(
                        intensity=intensity(eic_ms1[x,1]),
                        rt=rtime(eic_ms1[x,1]), mz=as.character(
                            rep(paste(round(ms1_peaks[x,"mz"],3),"m/z")))))
    df1 <- do.call("rbind", df1)
    
    p1 <- ggplot2::ggplot(df1[!is.na(df1$intensity),],
                            ggplot2::aes(x=rt,
                                        y=intensity, colour=as.factor(mz))) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::labs(x = "RT (s)", y = "Intensity (a.u.)",
                        colour = paste("Correlation >", cthres1)) +
        ggplot2::ggtitle(paste("Correlated EICs:", round(fmz,3),
                                "m/z",round(frt),"s")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5))
    
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
        ggplot2::ggtitle(paste("ISF Pseudo-MS/MS:", round(fmz,3),
                                "m/z",round(frt),"s")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5)) +
        ggplot2::ylim(0, max(df2[,2]) + 0.1*max(df2[,2])) +
        ggplot2::xlim(min(df2[,1])-50, max(df2[,1])+50) +
        ggplot2::labs(x = "m/z", y = "Intensity (a.u.)")
    
    if(length(aif) > 12){
        aif_ions <- match(aif[,"mz"], ms2_peaks[,"mz"])
    } else aif_ions <- match(aif["mz"], ms2_peaks[,"mz"])
    
    df3 <- lapply(aif_ions,
                    function(x) data.frame(
                                        intensity=intensity(eic_aif[x,1]),
                                        rt=xcms::rtime(eic_aif[x,1]), 
                                        mz=as.character(rep(
                                                paste(
                                                round(ms2_peaks[x,"mz"],3),
                                                "m/z")
                                                )
                                            )
                                        )
                )
    df3 <- do.call("rbind", df3)
    
    p3 <- ggplot2::ggplot(df3[!is.na(df3$intensity),],
        ggplot2::aes(x=rt, y=intensity, colour=mz)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::labs(x="RT (s)", y="Intensity (a.u.)",
                        colour=paste("Correlation >", cthres2)) +
        ggplot2::ggtitle(paste("Fragments EICs:", round(fmz,3),
                                "m/z",round(frt),"s")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5))
    
    if(length(aif) > 12){
        df4 <- as.data.frame(aif[,c("mz", "into")])
    } else df4 <- data.frame(mz=aif["mz"], into=aif["into"])
    
    p4 <- ggplot2::ggplot(df4,
        ggplot2::aes(x=mz, y=into, label = round(mz, 3))) +
        ggplot2::geom_segment(ggplot2::aes(xend = mz, yend=0),
                                color="red", lwd=0.5) +
        ggplot2::geom_text(size=3, angle=45, hjust=0, vjust=0) +
        ggplot2::ggtitle(paste("AIF Pseudo-MS/MS:", round(fmz,3),
                                "m/z",round(frt),"s")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +
        ggplot2::ylim(0, max(df4[,2]) + 0.1*max(df4[,2])) +
        ggplot2::xlim(min(df4[,1])-50, max(df4[,1])+50) +
        ggplot2::labs(x="m/z", y="Intensity (a.u.)")
    
    pdf(paste(DirPath, ExpName, "_", "pseudoMSMS_AIF_", round(fmz,3),
        "mz_", round(frt), "s", ".pdf", sep=""), height=8, width=12)
    gridExtra::grid.arrange(p1, p2, p3, p4, nrow=2)
    dev.off()
}

## plot EICs and insource fragmentation pseudo-MSMS using ggplot2 and gridExtra
plotISF <- function(fmz, frt, insource, ms1_peaks, 
                    eic_ms1, cthres1, DirPath, ExpName){
    if(length(insource) > 12){
        is_ions <- match(insource[,"mz"], ms1_peaks[,"mz"])
    } else is_ions <- match(insource["mz"], ms1_peaks[,"mz"])
    
    df1 <- lapply(is_ions,
                    function(x) data.frame(
                        intensity=intensity(eic_ms1[x,1]),
                        rt=xcms::rtime(eic_ms1[x,1]),
                        mz=as.character(
                                rep(paste(round(ms1_peaks[x,"mz"],3),"m/z")))))
    df1 <- do.call("rbind", df1)
    
    p1 <- ggplot2::ggplot(df1[!is.na(df1$intensity),],
                                ggplot2::aes(x=rt, y=intensity,
                                        colour=as.factor(mz))) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::labs(x="RT (s)", y="Intensity (a.u.)",
                        colour=paste("Correlation >", cthres1)) +
        ggplot2::ggtitle(paste("Correlated EICs:",
                                round(fmz,3),"m/z",round(frt),"s")) +
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
    
    if (length(insource) > 12){
        df2 <- as.data.frame(insource[,c("mz", "into")])
    } else df2 <- data.frame(mz=insource["mz"], into=insource["into"])
    
    p2 <- ggplot2::ggplot(df2,
                            ggplot2::aes(x=mz, y=into, label=round(mz,3))) +
        ggplot2::geom_segment(ggplot2::aes(xend = mz, yend=0),
                                color="red", lwd=0.5) +
        ggplot2::geom_text(size=3, angle=45, hjust=0, vjust=0) +
        ggplot2::theme_minimal() +
        ggplot2::ggtitle(paste("In-source MS spectrum:", round(fmz,3),"m/z", 
                                round(frt),"s")) +
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +
        ggplot2::ylim(0, max(df2[,2]) + 0.1*max(df2[,2])) +
        ggplot2::xlim(min(df2[,1])-50, max(df2[,1])+50) +
        ggplot2::labs(x="m/z", y="Intensity (a.u.)")
    
    pdf(paste(DirPath, ExpName, "_", "in-source_", round(fmz, 3),
                "mz_",round(frt),"s",".pdf", sep=""), height=8, width=12)
    gridExtra::grid.arrange(p1, p2, nrow=1)
    dev.off()
}

## save AIF or ISF pseudo-MS/MS as .csv and .mgf-------------------------------
savePseudoMSMS <- function(fmz, frt, xcmsF1, 
                            pseudo, specType, DirPath, ExpName){
    if(length(pseudo) > 12){
        df <- as.data.frame(pseudo[,c("mz", "into")])
    } else df <- data.frame(mz=pseudo["mz"], into=pseudo["into"])
    # save aif pseudo-MSMS to csv
    write.csv(df, paste(DirPath, ExpName, "_", "pseudoMSMS_", specType, "_", 
                        round(fmz, 3), "mz_",round(frt),"s", 
                        ".csv", sep=""))
    # save aif pseudo-MSMS in .mgf format
    # based on the definition from... 
    # https://fiehnlab.ucdavis.edu/projects/lipidblast/mgf-files
    # and
    # http://www.matrixscience.com/help/data_file_help.html
    fname <- paste(DirPath, ExpName, "_", "pseudoMSMS_", specType, "_", 
                    round(fmz, 3), "mz_",round(frt),"s", 
                    ".mgf", sep="")
    # get charge from polarity
    charge <- polarity(xcmsF1)
    charge <- unique(charge)
    if(charge == 1){
        charge <- "1+"
    } else if(charge == 0) charge <- "1-"
    cat("BEGIN IONS",paste("PEPMASS=", fmz, sep = ""),
        paste("CHARGE=", charge, sep = ""),
        paste("TITLE=", "Pseudo-MS/MS from", specType,"at", frt, "seconds"),
        sep = "\n", file = fname)
    for(i in seq_along(df$mz)){
        cat(paste(df$mz[i], " ", df$into[i], sep=""), "\n", file=fname, 
            append=TRUE)
    }
    cat("END IONS", "\n", file = fname, append=TRUE)
}