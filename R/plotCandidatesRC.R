#' Visual results of annotations based on LC-MS AIF processed using RAMClusteR.
#'
#' Plots pseudo-MS/MS spectra of matched fragments from RC objects.
#' Runs inside annotateRC function.
#'
#' @author Goncalo Graca & Yuheng (Rene) Cai (Imperial College London)
#'
#' @param fmz The m/z for the feature of interest.
#' @param frt Retention time in seconds for the feature of interest.
#' @param highCESpec MS2 peaks at the RT window of the feature of interest.
#' @param DatasetName Name of the data set processed using RAMClustR.
#' @param rankedCandidates List containing the ranked candidate annotations.
#' @param candidate Number of candidate annotations to plot.
#' @param DirPath Path to the folder where the plots will be saved.
#' @return Saves the plots of the candidate annotations as one pdf file:
#' a plot pseudo-MS/MS spectrum for the matched ions.
#' @export
plotCandidatesRC <- function(fmz, frt, highCESpec, DatasetName,
                             rankedCandidates, candidate, DirPath){

# get relevant information from the rankedCandidates object
metabolite <- rankedCandidates$rankedResult[candidate,"metabolite"]
ionType <- rankedCandidates$rankedResult[candidate,"ion.type"]
score <- rankedCandidates$rankedResult[candidate,"score"]
adductName <- paste(metabolite,ionType)
candidateMZ <- round(rankedCandidates$rankedResult[candidate,"mz.metabolite"],3)
MZerror <- round(rankedCandidates$rankedResult[candidate,"mz.error"],1)
specMatch <- rankedCandidates$rankedSpecMatch[[candidate]]
rnk <-  rankedCandidates$rankedResult[candidate,"rank"]

# plotting part
if(is.null(nrow(specMatch))) {
	  df1 <- as.data.frame(specMatch[c("mz", "into")])
} else df1 <- as.data.frame(specMatch[,c("mz", "into")])

p1 <- ggplot2::ggplot(df1,
                      ggplot2::aes(x = mz, y = into, label = round(mz, 3))) +
  ggplot2::geom_segment(ggplot2::aes(xend = mz, yend=0), color="red", lwd=0.5) +
  ggplot2::geom_text(size=3, angle=45, hjust=0, vjust=0) +
  ggplot2::ggtitle(paste("Feature ", round(fmz,4), "m/z", "_", round(frt,4),
      "s ", ", Rank ", rnk, " result: ", adductName, ", \u0394ppm =", MZerror,
      " ppm, ", "score = ", round(score, 2), sep = "")) +
  ggplot2::theme_minimal() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
  ggplot2::ylim(0, max(df1$into) + 0.1*max(df1$into)) +
  ggplot2::xlim(min(df1$mz) - 50, max(df1$mz) + 50) +
  ggplot2::labs(x = "m/z", y = "Intensity (a.u.)")

pdf(file = paste(DirPath, DatasetName, "_", round(fmz,3), "mz_", round(frt,3),
                 "_candidate_", candidate,".pdf", sep=""),
    width = 10, height = 5)
  gridExtra::grid.arrange(p1, nrow = 1)
  dev.off()
}
