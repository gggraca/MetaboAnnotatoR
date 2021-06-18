#' Visual results of single LC-MS AIF chromatogram-based annotations.
#'
#' Plot pseudo-MS/MS spectra of matched fragments and corresponding EICs.
#' Runs inside annotateAIF function.
#'
#' @author Goncalo Graca (Imperial College London)
#'
#' @param fmz The m/z for the feature of interest.
#' @param frt Retention time in seconds for the feature of interest.
#' @param iso Isotope "tag" to add to the results.
#' @param highCESpec MS2 peaks at the RT window of the feature of interest.
#' @param ms2_eic Object containing the EICs for the AIF features in the RT
#' window of the feature of interest;
#' @param SpName Sample name label.
#' @param rankedCandidates List containing the ranked candidate annotations.
#' @param candidate Number of candidate annotations to plot.
#' @param DirPath Path to the folder where the plots will be saved.
#' @return Saves the plots of the candidate annotations as one pdf file:
#' a plot of EICs and pseudo-MS/MS spectrum for the matched ions.
#' @export
plotCandidatesAIF <- function(fmz,
                              frt,
                              highCESpec,
                              ms2eic,
                              SpName,
                              rankedCandidates,
                              candidate,
                              DirPath){
# get relevant information from the rankedCandidates object
metabolite <- rankedCandidates$rankedResult[candidate,"metabolite"]
ionType <- rankedCandidates$rankedResult[candidate,"ion.type"]
score <- rankedCandidates$rankedResult[candidate,"score"]
adductName <- paste(metabolite,ionType)
candidateMZ <- round(rankedCandidates$rankedResult[candidate,"mz.metabolite"],3)
MZerror <- round(rankedCandidates$rankedResult[candidate,"mz.error"],1)
specMatch <- rankedCandidates$rankedSpecMatch[[candidate]]
rnk <- rankedCandidates$rankedResult[candidate,"rank"]

if(nrow(specMatch) >= 1){
# plotting part
  mz_idx <- match(specMatch[,1], highCESpec[,"mz"])
  df1 <- lapply(mz_idx, function(x) data.frame(
    intensity = xcms::intensity(ms2eic[[x]][1,1]),
    rt = xcms::rtime(ms2eic[[x]][1,1]),
    mz = as.character(rep(paste(round(highCESpec[x,"mz"],3),"m/z")))))
    df1 <- do.call("rbind", df1)

  p1 <- ggplot2::ggplot(df1[!is.na(df1$intensity),],
                        ggplot2::aes(x = rt, y = intensity, colour = mz)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(x = "RT (s)", y = "Intensity (a.u.)", colour = "fragments") +
    ggplot2::ggtitle(paste("Feature",round(fmz,3),"m/z",round(frt),"s,",
                    "Rank", rnk, "result:", adductName, ", \u0394ppm =",
                    MZerror, ", score =", round(score, 2))) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  df2 <- as.data.frame(specMatch[,c("mz", "into")])

  p2 <- ggplot2::ggplot(df2,
                        ggplot2::aes(x = mz, y = into, label = round(mz, 3))) +
    ggplot2::geom_segment(ggplot2::aes(xend = mz, yend=0),
                          color="red", lwd=0.5) +
    ggplot2::geom_text(size=3, angle=45, hjust=0, vjust=0) +
    ggplot2::ggtitle(paste("ions matched to", adductName)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::ylim(0, max(df2[,2]) + 0.1*max(df2[,2])) +
    ggplot2::xlim(min(df2[,1])-50, max(df2[,1])+50) +
    ggplot2::labs(x = "m/z", y = "Intensity (a.u.)")

  pdf(file = paste(DirPath,SpName,"_",round(fmz,3),"mz_",
                   round(frt,3),"_candidate_", candidate, ".pdf", sep=""),
      width = 10, height = 10)
    gridExtra::grid.arrange(p1, p2, nrow = 2)
  dev.off()
	} else NULL
}
