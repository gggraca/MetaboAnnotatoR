## Function plotCandidatesAIF to plot spectra of matched fragments and EICs
## Run only after running annotateAIF function
## Requires xcms version > 3.0 and .mzML raw data for EIC extraction
## Goncalo Graca 3/10/2020
## Adapted from plotCandidatesRaw to plot EIC objects from xcms version > 3.0
## 4/02/2019: change in plot title information
## 5/06/2019: adapted to work with compFrag2 and rankScore2 functions
## 24/09/2019: corrections to incorporate new names from rankedCandidates
plotCandidatesAIF <- function(fmz, frt, highCESpec, ms2eic, SpName, rankedCandidates, 
                              candidate, DirPath){
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
    df1 <- lapply(mz_idx, function(x) data.frame(intensity = intensity(ms2eic[[x]][1,1]), 
                                                        rt = rtime(ms2eic[[x]][1,1]), 
                                                        mz = as.character(rep(paste(round(highCESpec[x,"mz"],3),"m/z")))))
    df1 <- do.call("rbind", df1)
    
    p1 <- ggplot(df1[!is.na(df1$intensity),], aes(x = rt, y = intensity, colour = mz)) + geom_point() + geom_line() +
      labs(x = "RT (s)", y = "Intensity (a.u.)", colour = "fragments") +
      ggtitle(paste("Feature",round(fmz,3),"m/z",round(frt),"s,", "Rank", rnk, "result:", 
	  adductName, ", delta =", MZerror, ", score =", round(score, 2))) +
      theme(plot.title = element_text(hjust = 0.5))
    
    df2 <- as.data.frame(specMatch[,c("mz", "into")])
    
    p2 <- ggplot(df2, aes(x = mz, y = into, label = round(mz, 3))) + geom_segment(aes(xend = mz, yend=0), color="red", lwd=0.5) +
      geom_text(size=3, angle=45, hjust=0, vjust=0) + ggtitle(paste("ions matched to", adductName)) +
      theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) + ylim(0, max(df2[,2]) + 0.1*max(df2[,2])) + 
      xlim(min(df2[,1])-50, max(df2[,1])+50) + labs(x = "m/z", y = "Intensity (a.u.)")
    
    pdf(file = paste(DirPath,SpName,'_',round(fmz,3),"mz_",round(frt,3),'_candidate_',candidate,'.pdf', 
                     sep=''),width = 10, height = 10)
    grid.arrange(p1, p2, nrow = 2)
    dev.off()
	} else NULL
}