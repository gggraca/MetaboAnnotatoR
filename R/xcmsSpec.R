## Function xcmsSpec() ----------
## gets the low-colision-energy (MS1) spectrum or high-colision-energy spectrum (MS2) for any feature from an XCMS output
## fmz: feature m/z
## frt: feature RT in seconds
## xcmsObject: variable containing the XCMS processing object
## mztol: absolute tolerance for feature m/z search
## rttol: absolute tolerance for feature RT search
## highCE: high colision-energy option; if FALSE the "in-source" spectrum is returned
## requires XCMS package
## Goncalo Graca, g.gomes-da-graca@imperial.ac.uk

xcmsSpec <- function(fmz, frt, xcmsObject, rttol = 5, mztol = 0.01, highCE = TRUE){
  selection1 <- peaks(xcmsObject)[which(peaks(xcmsObject)[,"mz"] > (fmz - mztol) & 
                                         peaks(xcmsObject)[,"mz"] < (fmz + mztol) & 
                                         peaks(xcmsObject)[,"rt"] > (frt - rttol) & 
                                         peaks(xcmsObject)[,"rt"] < (frt + rttol)),]
  
  if(is.null(dim(selection1))) selection1 <- rbind(selection1, c(0,0))
  
  # get index of the sample with highest feature intensity
  # index corresponding to the highCE scans for the sample with highest feature intensity:
  if(highCE) idx <- as.numeric(selection1[which.max(selection1[,"into"]), "sample"] + 1)
  
  # index corresponding to the lowCE scans for the sample with highest feature intensity:
  if(!highCE) idx <- as.numeric(selection1[which.max(selection1[,"into"]), "sample"])
  
  # get all peaks from the selected sample around the feature rt
  selection2 <- which(peaks(xcmsObject)[,"sample"] == idx & 
                      peaks(xcmsObject)[,"rt"] > (frt - rttol) & 
                        peaks(xcmsObject)[,"rt"] < (frt + rttol))
  peakSelection <- peaks(xcmsObject)[selection2,]
  
  # check if peakSelection is vector or data frame and store result
  if(is.null(dim(peakSelection))) {
    result <- peakSelection[c("mz","into")]
    } else result <- peakSelection[, c("mz","into")]
  
  return(result)
}