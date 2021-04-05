#' Extracts the low-collision-energy (MS1) or high-collision-energy (MS2)
#' features (pseudo-spectra) for any feature from an XCMS output object.
#'
#' Given a feature of interest (m/z RT pair), the function will extract the
#' low-collision-energy (MS1) or high-collision-energy features (MS2)
#' at the same RT window of the feature of interest from an XCMS output object
#' obtained by processing both functions together from a set of LC-MS
#' chromatograms.
#'
#' @author Goncalo Graca & Yuheng (Rene) Cai (Imperial College London)
#'
#' @param fmz The m/z for the feature of interest.
#' @param frt Retention time in seconds for the feature of interest.
#' @param xcmsObject Variable containing the XCMS processing object.
#' @param mztol Absolute tolerance for feature m/z search in Da.
#' @param rttol Tolerance for feature RT search in seconds.
#' The default (5 s) only applies to UPLC/UHPLC data.
#' @param highCE Logic value. If TRUE the high collision-energy is extracted,
#' otherwise if FALSE the "in-source" spectrum is returned.
#' @returns A data frame with ions (m/z and intensity) from the
#' high collision-energy or low collision-energy features found at the same
#' RT window as the feature of interest.
#' @export
xcmsSpec <- function(fmz,
                     frt,
                     xcmsObject,
                     mztol = 0.01,
                     rttol = 5,
                     highCE = TRUE){

  selection1 <- peaks(xcmsObject)[which(
    xcms::peaks(xcmsObject)[,"mz"] > (fmz - mztol) &
    xcms::peaks(xcmsObject)[,"mz"] < (fmz + mztol) &
    xcms::peaks(xcmsObject)[,"rt"] > (frt - rttol) &
    xcms::peaks(xcmsObject)[,"rt"] < (frt + rttol)
    ),]

  if(is.null(dim(selection1))) selection1 <- rbind(selection1, c(0,0))

  # get index of the sample with highest feature intensity
  # index corresponding to the highCE scans for the sample
  # with highest feature intensity:
  if(highCE) {
    idx <- as.numeric(selection1[which.max(selection1[,"into"]), "sample"] + 1)
  }
  # index corresponding to the lowCE scans for the sample with
  # highest feature intensity:
  if(!highCE) {
    idx <- as.numeric(selection1[which.max(selection1[,"into"]), "sample"])
  }
  # get all peaks from the selected sample around the feature rt
  selection2 <- which(xcms::peaks(xcmsObject)[,"sample"] == idx &
                      xcms::peaks(xcmsObject)[,"rt"] > (frt - rttol) &
                      xcms::peaks(xcmsObject)[,"rt"] < (frt + rttol))
  peakSelection <- xcms::peaks(xcmsObject)[selection2,]

  # check if peakSelection is vector or data frame and store result
  if(is.null(dim(peakSelection))) {
    result <- peakSelection[c("mz","into")]
    } else result <- peakSelection[, c("mz","into")]

  return(result)
}
