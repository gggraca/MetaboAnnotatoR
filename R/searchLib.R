#' Searches candidate metabolites from the fragments libraries by m/z
#'
#' Search candidate metabolites from the fragments libraries
#' for a given feature using m/z and RT (if metabolite RTs are known).
#' If no match is found in the "parent" ions, for instance in the case of a
#' feature corresponding to an in-source fragment, fragments are then searched.
#' Fragment ions will only be considered if the parent ion is present in the
#' same pseudo-MS spectrum (MS1).
#'
#' @author Goncalo Graca & Yuheng (Rene) Cai (Imperial College London)
#'
#' @param libraries List object containing all loaded libraries.
#' @param libfiles Path to the libraries files.
#' @param fmz The m/z for the feature of interest.
#' @param frt Retention time in seconds for the feature of interest.
#' @param tolerance Tolerance for m/z candidate search in ppm.
#' @param RTs Optional metabolites classes Retention Times in seconds.
#' Default value is "none".
#' @param inSourceSpec Data frame containing the pseudo-MS spectrum (MS1).
#' This will be used to check for the "parent" ion when the feature of interest
#' is matched to a fragment (in-source fragment).
#' @return A list of data frames containing the candidates from the fragment
#' libraries which will be used in the pseudo-MS/MS to fragment matching step.
#' @export
searchLib <- function(libraries,
                      libfiles,
                      fmz,
                      frt,
                      tolerance = 25,
                      RTs,
                      inSourceSpec){

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
      # extract candidates; row with scores added again:
      candidates[[j]] <- rbind(tlib[which(MZerr <= tolerance),],scores)
      # add libraries index for fragment comparison:
      names(candidates)[j] <- as.character(i)
    } else {
	# check for fragments
		MZerr <- abs(tlib[,3:dim(tlib)[2]] - fmz) * 1e6 / fmz # ppm
			if (sum(MZerr <= tolerance) > 0){
				# check if the parent ion is present in the in source ions list
			  # store candidates in a temporary object:
				tmp <- tlib[which(MZerr < tolerance, arr.ind = TRUE)[,1],]
				# creates a vector of indexed entries that will be discarded (if any)
				idx <- seq(1:nrow(tmp))
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
				# row with scores added again
				candidates[[j]] <- rbind(tmp,scores)
				# add libraries index for fragment comparison
				names(candidates)[j] <- as.character(i)
			} else next
		} else next
	}
	}, error=function(e){NULL}) # removed: {message("",conditionMessage(e))})
	}
return(candidates)
}
