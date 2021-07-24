#' Function to generate metabolite database entries from MS/MS spectra 
#' obtained from from public databases, stored as a .txt object containing
#' m/z and intensity values, and read imported into R as matrix
#' 
#' @author Goncalo Graca (Imperial College London)
#' 
#' @param specObject a matrix or data frame object containing the MS/MS spectrum 
#' arranged in two columns: 'mz' and  'intensity'. Intensity can be given in 
#' absolute or relative scale.
#' @param name Metabolite name.
#' @param adduct Type of adduct fragmented.
#' @param tmz Adduct m/z.
#' @param noise Noise intensity threshold expressed as a ratio to the peak with
#' the highest intensity.
#' @param mpeaksScore The occurrence score to be attributed to the most intense 
#' peaks of the MS/MS spectrum which should correspond to the most characteristic
#' fragmentation ions from the metabolite (or 'marker' peaks). These will be the 
#' peaks above 'mpeaksThres' value. This score is divided by the number of peaks 
#' above 'mpeaksThres' threshold. By default this value is defined at 0.9, 
#' which means that peaks below 'mpeaksThres' threshold will be given an 
#' occurrence score of 0.1, so that the sum of all fragment occurrence scores is 1.
#' @param mpeaksThres Intensity threshold to select peaks of the MS/MS spectrum 
#' considered to be highest intensity, expressed as a ratio to the peak 
#' with the highest intensity.
#' @return A .csv file containing fragment and parent m/z values and corresponding 
#' occurrence scores.
#' @export
genFragEntry <- function(specObject,
                         name,
                         adduct,
                         tmz,
                         filename,
                         noise = 0.005,
                         mpeaksScore = 0.9, 
                         mpeaksThres = 0.1,
                         mzTol = 0.01) {
# sort and filter m/z values find maximum intensity peak for spectrum 
# normalisation and normalise spectrum
  specObject <- specObject[order(-specObject[,1]),]
  specObject <- specObject[which(specObject[,1] <= tmz+mzTol),]
	norm.specObject <- specObject
	norm.specObject[,2] <- specObject[,2]/specObject[which.max(specObject[,2]),2]
# denoise spectrum
	denoised.spec <- norm.specObject[which(norm.specObject[,2] > noise),]
# locate if parent m/z is present in the list
	p <- which.min(abs(denoised.spec[,1]-tmz))
	if(length(p)==0) denoised.spec <- rbind(c(tmz,0),denoised.spec)
	if(length(p)==1) denoised.spec[p,1] <- tmz
# define marker peaks and attribute score
	scores <- rep(0,nrow(denoised.spec))
	idx <- which(denoised.spec[,2] >= mpeaksThres)
	scores[idx] <- mpeaksScore/length(idx)
# attribute scores to the remaining peaks
	idx2 <- which(denoised.spec[,2] < mpeaksThres)
	scores[idx2] <- (1-mpeaksScore)/length(idx2)
# check if score adds to 1, if not recalculate scores
	if(sum(scores) < 1){
	  scores[idx] <- 1/length(idx)
	} else NULL
# save entry as .csv
	result <- rbind(denoised.spec[,1],scores)
	frag <- rep(NA,(length(scores) - 1))
	for (i in 1:length(scores)-1){
		frag[i] <- paste('fragment',i,sep='')
	}
	colnames(result) <- c(adduct, frag)
	rownames(result) <- c(name,'scores')
	write.csv(result, paste(filename,'.csv', sep = ''), row.names = TRUE)
}