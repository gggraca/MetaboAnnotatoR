#' Function to generate metabolite database entries from MS/MS spectra 
#' obtained from from an .msp file
#' 
#' @author Goncalo Graca (Imperial College London)
#' 
#' @param msp_file an MS/MS spectral library for spectra from one or both polarities
#' @param library_name Custom library name under which POS and NEG folders will be created 
#' and where the respective libray entries will be stored 
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
#' @param mzTol Absolute tolerance for m/z search in Da.
#' @return A .csv file containing fragment and parent m/z values and corresponding 
#' occurrence scores.
#' @examples 
#' mspToLib("MassBank_example.msp")
#' @export
mspToLib <- function(msp_file,
                     library_name = "Custom",
                     noise = 0.005,
                     mpeaksScore = 0.9, 
                     mpeaksThres = 0.1) {
  
  # create folder to store library
  if(dir.exists("./Libraries")){
    dir.create(paste("./Libraries/",library_name, sep=""), showWarnings = FALSE)
    dirPath <- paste("./Libraries/",library_name, sep="")
  } else {
    dir.create("./Libraries/", showWarnings = FALSE)
    dir.create(paste("./Libraries/",library_name, sep=""), showWarnings = FALSE)
    dirPath <- paste("./Libraries/",library_name, sep="")
  }
    
  # read msp file
  m <- readLines(msp_file, warn = FALSE)
  
  #get names of all metabolites
  n <- grep("Name:", m)
  names <- unlist(lapply(n, function(x) substring(m[x], 7, nchar(m[x]))))
  
  # get number of peaks per MS/MS record
  np <- grep("Num Peaks:", m)
  npeaks <- unlist(lapply(np, function(x) substring(m[x], 11, nchar(m[x]))))
  npeaks <- as.numeric(npeaks)
  
  # get all ion modes
  im <- grep("Ion_mode:", m)
  ion_modes <- unlist(lapply(im, function(x) substring(m[x], 11, nchar(m[x]))))
  
  getMSPdetails <- function(x){
    # get block of text
    c <- m[n[x]:np[x]]
    
    # get precursor m/z
    pmz <- grep("PrecursorMZ:", c)
    if(length(pmz) == 0) precursor_mz <- NA else {
      precursor_mz <- substring(c[pmz], 14, nchar(c[pmz]))
      precursor_mz <- as.numeric(precursor_mz)
    }
    
    # get ion mode
    im <- grep("Ion_mode:", c)
    if(length(im) == 0) ion_mode <- NA else {
      ion_mode <- substring(c[im], 11, nchar(c[im]))
    }
    
    # get Precursor_type:
    pt <- grep("Precursor_type:", c)
    if(length(pt) == 0) ptype <- NA else {
      ptype <- substring(c[pt], 17, nchar(c[pt]))
    }
    
    # get MS/MS spectrum
    s <- NULL
    for(i in 1:npeaks[x]){
      tmp <- strsplit(m[np[x]+i], split = " ")
      s <- c(s, as.numeric(tmp[[1]]))
    }
    spec <- matrix(s, ncol = 2, byrow = TRUE)
    
    result <- list(metabolite = names[x],
                   precursor = precursor_mz,
                   type = ptype,
                   ion_mode = ion_mode,
                   MSMS = spec)
    return(result)
  }
  
  libs <- lapply(1:length(n), function(x) getMSPdetails(x))
  
  for(i in 1:length(libs)){
    name <- libs[[i]]$metabolite
    ion_mode <- libs[[i]]$ion_mode
    adduct <- libs[[i]]$type
    tmz <- libs[[i]]$precursor
    filename <- paste(name,".csv", sep = "")
    specObject <- libs[[i]]$MSMS
    
    # sort and filter m/z values find maximum intensity peak for spectrum
    if(nrow(specObject)>1){
      specObject <- specObject[order(-specObject[,1]),]
    } else NULL
    
    # normalise spectrum
    norm.specObject <- specObject
    norm.specObject[,2] <- specObject[,2]/specObject[which.max(specObject[,2]),2]
    
    # denoise spectrum
    if(nrow(specObject)>1){
      denoised.spec <- norm.specObject[which(norm.specObject[,2] > noise),]
    } else denoised.spec <- norm.specObject
    
    if(is.vector(denoised.spec)) denoised.spec <- matrix(denoised.spec, nrow = 1)
      
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
    if(nrow(specObject)>1){
      result <- rbind(denoised.spec[,1],scores)
      frag <- rep(NA,(length(scores) - 1))
      for (i in 1:length(scores)-1){
        frag[i] <- paste('fragment',i,sep='')
      }
      colnames(result) <- c(adduct, frag)
      rownames(result) <- c(name,'scores')
    } else if(nrow(specObject) == 1) {
      result <- rbind(specObject[,1],1)
      colnames(result) <- adduct
      rownames(result) <- c(name,'scores')
    } else NULL
    if(ion_mode=="POSITIVE") {
      dir.create(paste(dirPath,"/POS/",sep=""), showWarnings = FALSE)
      targetPath <- paste(dirPath,"/POS/", filename, sep = "")
    } else if(ion_mode=="NEGATIVE"){
      dir.create(paste(dirPath,"/NEG/",sep=""), showWarnings = FALSE)
      targetPath <- paste(dirPath,"/NEG/", filename, sep = "")
    }
    write.csv(result, targetPath, row.names = TRUE)
  }
}
