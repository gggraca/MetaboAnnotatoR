## Function compFrag() ----------
### Goncalo Graca & Yuheng Cai (Imperial College London)
## g.gomes-da-graca@imperial.ac.uk
# compare high-collision-energy spectrum m/z with fragments databases
# fragments are selected from accurate MS search and expected retention time
# reduced mz.tolerance from compFrag to 0.01

compFrag <- function(candidate, lib, fmz, frt, iso, highCESpec, pseudoSpec, tolerance = 25, maxMZdiff = 0.01, 
                     matchWeight = 0.5, NoMatchWeight = 0.5, useMZerrorWeight = TRUE, additional = TRUE){
  
  # check if highCESpec is vector...important when working with RC objects
  # if true stack extra row
  if(is.null(dim(highCESpec))){
    highCESpec <- rbind(highCESpec, c(0,0))
  }
  
  # calculate MZerrorWeight 
  MZerrorWeight <- 1 - matchWeight
  if(!useMZerrorWeight) matchWeight <- 1
  
  # strip the scores from the last row of the candidate data frame
  scores <- candidate[nrow(candidate),]
  candidate <- candidate[-nrow(candidate),]
  
  # temporary results table setup
  tempResult<-data.frame(cbind(rep(fmz,dim(candidate)[1]),rep(frt,dim(candidate)[1])))
  colnames(tempResult)[1:2] <- c("feature.mz","feature.rt")
  tempResult[c("metabolite","feature.type","ion.type","isotope","mz.metabolite","matched.mz",
  "mz.error","pseudoMSMS","fraction","score")] <- NA
  
  # type of ion isotope
  if (iso==0) {
	tempResult$isotope <- "M+0"
  } else if (iso==1) {
	tempResult$isotope <- "M+1"
  } else if (iso==2) {
	tempResult$isotope <- "M+2"
  } else if (iso==3) {
	tempResult$isotope <- "M+3"
  }

  # candidate name and m/z
  tempResult$metabolite <- candidate[,1]
  tempResult$mz.metabolite <- candidate[,2]
  
  # pseudoMSMS flag
  if (is.null(pseudoSpec) | length(pseudoSpec) == 0) {
	    tempResult$pseudoMSMS <- 'FALSE'
	  } else { tempResult$pseudoMSMS <- 'TRUE'
  }
  
  # spec match and final results table setup 
  specMatch <- list()
  result <- NULL
  
  # feature type by comparison with candidate
  for (i in 1:dim(candidate)[1]){
  # 4 is the the smallest difference between parent and fragment assuming isotopes up to M+3:
	if (abs(fmz-tempResult$mz.metabolite[i])<= 4) { tempResult$feature.type <- "parent"
		} else tempResult$feature.type <- "fragment"
  }
  
  for (i in 1:dim(candidate)[1]){
	# m/z error of feature
	correctedMZ <- fmz - iso*1.0034
    tempResult$mz.error[i] <- min(abs(candidate[i,2:dim(candidate)[2]] - correctedMZ) * 1e6 / correctedMZ)
	# m/z matched parent or fragment
	tempResult$matched.mz[i] <- candidate[i, which.min(abs(candidate[i,2:dim(candidate)[2]] - correctedMZ) * 1e6 / correctedMZ) + 1]
	tempResult$ion.type[i] <- colnames(candidate)[which.min(abs(candidate[i,2:dim(candidate)[2]] - correctedMZ) * 1e6 / correctedMZ) + 1]
	
	# score and matches counters
    score <- 0
    count <- 0
	nonMatched <- 0
  
  # store matched m/z and corresponding intensity and also non-matched mz for second search in hce spectrum
    mz <- vector()
    into <- vector()
	nonMatchedPseudo <- vector()
    nonMatchedScores <- vector()
	
    # matching (including parent ion m/z)
    for (j in 2:dim(candidate)[2]){
   	# comparison between Library entry and pseudoMSMS spectrum
      if (!is.null(pseudoSpec)) { # if (!is.null(pseudoSpec) | !is.na(pseudoSpec)){
        if (is.vector(pseudoSpec)) pseudoSpec <- as.data.frame(t(pseudoSpec)) else NULL # to avoid craches when pseudoSpec has only one entry
        if (any(abs(pseudoSpec[,'mz'] - candidate[i,j]) < maxMZdiff) | any(abs(pseudoSpec[,'mz']-candidate[i,j]) < maxMZdiff)){
          count <- count + 1
          score <- score + scores[,j]
		  mz[count] <- pseudoSpec[which.min(abs(pseudoSpec[,1] - candidate[i,j])), "mz"]
    into[count] <- pseudoSpec[which.min(abs(pseudoSpec[,1] - candidate[i,j])), "into"]
        } else {
		    nonMatched <- nonMatched + 1
		    nonMatchedPseudo[nonMatched] <- candidate[i,j]
		    nonMatchedScores[nonMatched] <- scores[,j]
		}
	# if pseudoMSMS spectrum does not exist, then compare with hce spectrum	
      } else if (is.null(pseudoSpec)) { # removed to avoid crashes with annotateRC : if (is.null(pseudoSpec) | length(pseudoSpec)==0 | is.na(pseudoSpec))
		    if (any(abs(highCESpec[,"mz"] - candidate[i,j]) < maxMZdiff)){
        count <- count + 1
        score <- score + scores[,j] * NoMatchWeight
        # store matched m/z and corresponding intensity
        mz[count] <- highCESpec[which.min(abs(highCESpec[,"mz"] - candidate[i,j])), "mz"]
        into[count] <- highCESpec[which.min(abs(highCESpec[,"mz"] - candidate[i,j])), "into"]
		    }
	    }
    }    
	# if some fragments are not matched to pseudoMSMS, try to match the with hce spectrum; consider giving less weigth to the matches here
	if(additional & !is.null(pseudoSpec)){ # !is.null(pseudoSpec) to avoid crashes with annotateRC
	  second.count <- 0
	if (nrow(pseudoSpec) > 0 & length(nonMatchedPseudo) > 0){
		for (k in 1:length(nonMatchedPseudo)){
			if (any(abs(highCESpec[,"mz"] - nonMatchedPseudo[k]) < maxMZdiff)){
			count <- count + 1
			score <- score + nonMatchedScores[k] * NoMatchWeight
			second.count <- second.count + 1 # add second count to denominator in score calculation
			# store additional matched m/z and corresponding intensity
			mz[count] <- highCESpec[which.min(abs(highCESpec[,"mz"] - nonMatchedPseudo[k])), "mz"]
			into[count] <- highCESpec[which.min(abs(highCESpec[,"mz"] - nonMatchedPseudo[k])), "into"]
		}
	 }
	}
	}
	# store fraction, score and matched spectrum
	#if (tempResult$mz.error[i]==0) {
	#tempResult$score[i] <- matchWeight*score	# to prevent division by zero
	#} else tempResult$score[i] <- matchWeight*score + MZerrorWeight*(1/tempResult$mz.error[i])
	if(score > 0) {
		if(useMZerrorWeight){
			if (tempResult$mz.error[i] < 1) {
			tempResult$score[i] <- MZerrorWeight + matchWeight*score	# to prevent division by zero
		} else tempResult$score[i] <- matchWeight*score + MZerrorWeight*(1/tempResult$mz.error[i])
		tempResult$fraction[i] <- paste(count," of ",(length(candidate[i,])-1))
		if(tempResult$score[i] > 0){
			specMatch[[i]] <- data.frame("mz" = mz, "into" = into)
			names(specMatch)[i] <- as.character(candidate[i,1]) # save name of candidate
			result <- rbind(result,tempResult[i,])
			}
		}
		if(!useMZerrorWeight){
			tempResult$score[i] <- matchWeight*score
			tempResult$fraction[i] <- paste(count," of ",(length(candidate[i,]) - 1))
			if(tempResult$score[i] > 0){
				specMatch[[i]] <- data.frame("mz" = mz, "into" = into)
				names(specMatch)[i] <- as.character(candidate[i,1]) # save name of candidate
			result <- rbind(result, tempResult[i,])
		} else next
	}  
  } else next
  }
 return(list("result" = result, "specMatch" = specMatch))
}