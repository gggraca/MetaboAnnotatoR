#' Function to generate metabolite database entries from MS/MS spectra 
#' obtained from from an .msp file
#' 
#' @author Goncalo Graca (Imperial College London)
#' 
#' @param msp_file an MS/MS spectral library for spectra from one or both polarities
#' @return A .csv file containing fragment and parent m/z values and corresponding 
#' occurrence scores.
#' @export
mspToLib <- function(msp_file) {
m <- readLines(msp_file)

#get names of all metabolites
n <- grep("Name:", m)
names <- unlist(lapply(n, function(x) substring(m[x], 7, nchar(m[x]))))

# get number of peaks per MS/MS record
np <- grep("Num Peaks:", m)
names <- unlist(lapply(n, function(x) substring(m[x], 7, nchar(m[x]))))
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
  # use different strategy...perhaps read line and split elements
  s <- NULL
  for(i in 1:npeaks[x]){
    tmp <- strsplit(m[np[x]+i], split = " ")
    s <- c(s, as.numeric(tmp[[1]]))
  }
  spec <- matrix(s, ncol = 2, byrow = TRUE)
  # spec <- scan("MassBank_NIST.msp", skip = np[x], nlines = npeaks[x], quiet = TRUE)
  
  result <- list(metabolite = names[x],
                 precursor = precursor_mz,
                 type = ptype,
                 ion_mode = ion_mode,
                 MSMS = spec)
  return(result)
}

libs <- lapply(1:length(n), function(x) getMSPdetails(x))

# all records to the same folder...
# needs changes sort records automatically...
l <- lapply(1:length(libs), function(x) genFragEntry(specObject = libs[[x]]$MSMS,
                                                     name = libs[[x]]$metabolite,
                                                     adduct = libs[[x]]$type,
                                                     tmz = libs[[x]]$precursor,
                                                     filename = paste(libs[[x]]$metabolite, ".csv", sep = ""),
                                                     noise = 0.005,
                                                     mpeaksScore = 0.9, 
                                                     mpeaksThres = 0.1,
                                                     mzTol = 0.01))
}
