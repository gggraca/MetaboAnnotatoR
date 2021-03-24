## Function searchLib() ----------
## Goncalo Graca & Yuheng Cai (Imperial College London)
## g.gomes-da-graca@imperial.ac.uk
# search matched metabolites for a given feature with rt, mz
# Replace first part of the function to read the same information from a .csv or .txt file
# If the file is a .csv file, care must be taken to account for the commas
# Add for cycle to go through all lines and verify correspondence between feature and lipid class RT
# 18/01/2019: tryCatch function added to bypass errors in fragment search 
# 29/01/2019: RT information stored in an external file read into the object 'RTs'
# 15/06/2019: extra coded added filters candidates not present in the retention time window of the feature of interest
# Important note: replace peaksF1 with full MS info to work with RC objects
# 5/11/2019: Function adapted to cope with the additional scores row in the libraries
# 23/01/2020: Instead of search for parents in the raw data retention time window, search is made on the in source ions only to reduce the number of false positives
# 28/01/2020: Corrections added to prevent duplication of candidate selection
# 28/01/2020: parent filtration based on in source fragments improved

searchLib <- function(libraries, libfiles, fmz, frt, tolerance = 25, RTs, inSourceSpec){
  # narrow down candidate classes by RT (min)
  classgroup <- NULL
  
  if(is.character(RTs)) ilib <- seq(1,length(libfiles),1) 
  if(!is.character(RTs)){
  for(i in 1:dim(RTs)[1]){
		if(frt/60>RTs[i,1] & frt/60<RTs[i,2]) classgroup <- RTs[i,3:dim(RTs)[2]]
		}
	classgroup <- classgroup[-which(classgroup == "")]
	ilib <- unlist(lapply(classgroup, grep, x = libfiles))
	ilib <- unique(ilib)
  }
# obtain candidates by mz
  candidates <- list()
  j <- 0 # index of candidates
  for (i in ilib){
    tryCatch({
	# strip scores row
	tlib <- libraries[[i]][-nrow(libraries[[i]]),]
	scores <- libraries[[i]][nrow(libraries[[i]]),]
	MZerr <- abs(tlib[,2] - fmz) * 1e6 / fmz # ppm ; row with scores removed
    if (sum(MZerr < tolerance) > 0){
      j <- j+1
      candidates[[j]] <- rbind(tlib[which(MZerr <= tolerance),],scores) # extract candidates; row with scores added again
      names(candidates)[j] <- as.character(i) # add libraries index for fragment comparison 
    } else {
	# check for fragments
		MZerr <- abs(tlib[,3:dim(tlib)[2]] - fmz) * 1e6 / fmz # ppm
			if (sum(MZerr <= tolerance) > 0){
				# check if the parent ion is present in the in source ions list
				tmp <- tlib[which(MZerr < tolerance, arr.ind = TRUE)[,1],] # store candidates in a temporary object
				idx <- seq(1:nrow(tmp)) # create a vector of indexed for the entries that will be descarted if any
				for (k in 1:nrow(tmp)){
					res <- which(abs(inSourceSpec[,'mz']-tmp[k,2])<0.01)
					if(length(res)==0) {
					idx[k] <- idx[k]*-1
					} else next
				}
			#tmp <- tmp[idx,]
			tmp <- tmp[which(idx >0),]
			if(nrow(tmp)>0) {
				j<-j+1
				candidates[[j]] <- rbind(tmp,scores) # row with scores added again
				names(candidates)[j] <- as.character(i) # add libraries index for fragment comparison
			} else next		
		} else next
	}
	}, error=function(e){NULL}) # removed: {message("",conditionMessage(e))}) 
	}	
return(candidates)
}
