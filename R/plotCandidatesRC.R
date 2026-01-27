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
#' @param candidate Index of the candidate in the rankedCandidates list to 
#' plot annotations from.
#' @param DirPath Path to the folder where the plots will be saved.
#' @return Saves the plots of the candidate annotations as one pdf file:
#' a plot pseudo-MS/MS spectrum for the matched ions.
#' @examples
#' setwd(tempdir())
#' # create a results data frame for 3 hypothetical metabolite candidates
#' # for feature 152.0723 m/z and 125 s
#' fmz <- 152.0723
#' frt <- 125
#' result <- data.frame(metabolite=c("Met A", "Met B", "Met C"), 
#' feature.type=rep("parent",3), ion.type=rep("[M+H]+"),
#' isotope=rep("M+0",3), mz.metabolite=rep(152.0723, 3), 
#' matched.mz=rep(152.0706, 3), mz.error=rep(11, 3),
#' pseudoMSMS=rep("TRUE", 3), fraction=c("2 of 5", "4 of 5","3 of 5"), 
#' score=c(0.4, 0.9, 0.6))
#' # create entries for spectra matches of the hypothetical metabolite candidates
#' specMatch <- list()
#' specMatch$`Met A` <- data.frame(mz=c(152.0716, 134.0611, 59.0489, 65.0389, 66.0427), 
#' into=c(432, 592, 2092, 4836, 832))
#' specMatch$`Met B` <- data.frame(mz=c(152.0716, 134.0611, 110.0622, 109.0523, 59.0489), 
#' into=c(65236, 4480, 30696, 448, 432)) 
#' specMatch$`Met C` <- data.frame(mz=c(152.0716, 134.0611, 110.0622, 109.0523, 93.0569), 
#' into=c(65236, 4480, 30696, 464, 804))
#' # create high collision energy spectrum as data frame
#' highCESpec <- data.frame(mz=c(59.0489, 65.0389, 66.0427, 67.0550, 70.0659, 73.0762, 
#' 82.0658, 92.0498, 93.0355, 93.0569, 109.0523, 110.0622, 111.0452, 111.0647, 112.0476, 
#' 121.0408, 134.0611, 136.0762, 152.0716, 154.0781), 
#' into=c(3228, 8696, 564, 1004, 432, 592, 2092, 4836, 832, 560, 448, 30696, 8516, 3400, 
#' 464, 804, 4480, 368, 65236, 464))
#' # run rankScore function and get the ranked result
#' rankedCandidates <- rankScore(result, specMatch)
#' # Plot the first candidate in the candidates list
#' plotCandidatesRC(fmz, frt, highCESpec, DatasetName="LCMS", rankedCandidates, 
#' candidate=1, DirPath=getwd())
#' @export
plotCandidatesRC <- function(fmz, frt, highCESpec, DatasetName,
                                rankedCandidates, candidate, DirPath){

    # get relevant information from the rankedCandidates object
    metabolite <- rankedCandidates$rankedResult[candidate,"metabolite"]
    ionType <- rankedCandidates$rankedResult[candidate,"ion.type"]
    score <- rankedCandidates$rankedResult[candidate,"score"]
    adductName <- paste(metabolite,ionType)
    candidateMZ <- round(
        rankedCandidates$rankedResult[candidate,"mz.metabolite"],3)
    MZerror <- round(rankedCandidates$rankedResult[candidate,"mz.error"],1)
    specMatch <- rankedCandidates$rankedSpecMatch[[candidate]]
    rnk <-  rankedCandidates$rankedResult[candidate,"rank"]

    # plotting part
    if(is.null(nrow(specMatch))) {
        df1 <- as.data.frame(specMatch[c("mz", "into")])
    } else df1 <- as.data.frame(specMatch[,c("mz", "into")])

    p1 <- ggplot2::ggplot(df1,
                            ggplot2::aes(x=mz, y=into, label=round(mz, 3))) +
        ggplot2::geom_segment(ggplot2::aes(xend=mz, yend=0), 
                                color="red", lwd=0.5) +
        ggplot2::geom_text(size=3, angle=45, hjust=0, vjust=0) +
        ggplot2::ggtitle(paste("Feature ", round(fmz,4), "m/z", "_", 
                                round(frt,4), "s ", ", Rank ", rnk, 
                                " result: ", adductName, ", mz.error = ", 
                                MZerror, " ppm, ", "score = ", 
                                round(score, 2), sep="")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +
        ggplot2::ylim(0, max(df1$into) + 0.1*max(df1$into)) +
        ggplot2::xlim(min(df1$mz)-50, max(df1$mz)+50) +
        ggplot2::labs(x="m/z", y="Intensity (a.u.)")

    pdf(file=paste(DirPath, DatasetName, "_", round(fmz,3), "mz_", 
                    round(frt,3), "_candidate_", candidate,".pdf", sep=""), 
        width = 10, height = 5)
    gridExtra::grid.arrange(p1, nrow = 1)
    dev.off()
}
