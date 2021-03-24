plotCandidatesRC <- function(fmz, frt, highCESpec, DatasetName, rankedCandidates, 
                              candidate, DirPath){
 
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
    
#	df2 <- as.data.frame(highCESpec[,c("mz", "into")])

      p1 <- ggplot(df1, aes(x = mz, y = into, label = round(mz, 3))) + 
      geom_segment(aes(xend = mz, yend=0), color="red", lwd=0.5) +
      geom_text(size=3, angle=45, hjust=0, vjust=0) + 
      ggtitle(paste("Feature ", round(fmz,4), "m/z", "_", round(frt,4), 
      "s ", ", Rank ", rnk, " result: ", adductName, ", delta =", MZerror, 
      " ppm, ", "score = ", round(score, 2), sep = "")) +
      theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) + 
      ylim(0, max(df1$into) + 0.1*max(df1$into)) + 
      xlim(min(df1$mz) - 50, max(df1$mz) + 50) + labs(x = "m/z", y = "Intensity (a.u.)")
    
    # p2 <- ggplot(df2, aes(x = mz, y = into, label = round(mz, 3))) + 
    #   geom_segment(aes(xend = mz, yend=0), color="red", lwd=0.5) +
    #   geom_text(size = 3, angle = 45, hjust = 0, vjust = 0) + 
    #   ggtitle(paste("AIF spectrum", round(frt, 4) , " +/- 5 s")) +
    #   theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) +
    #   ylim(0, max(df2$into) + 0.1*max(df2$into)) + 
    #   xlim(min(df2$mz), max(df2$mz)+50) + labs(x = "m/z", y = "Intensity (a.u.)")
    
    pdf(file = paste(DirPath,DatasetName,'_',round(fmz,3),"mz_",round(frt,3),'_candidate_',
                    candidate,'.pdf', sep=''), width = 10, height = 5)
    grid.arrange(p1, nrow = 1)
    dev.off()
}