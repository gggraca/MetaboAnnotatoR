#' Function checkIsotope
#' checks the type of isotope of the input feature
#' returns the "tag" of the isotope from an isotopic series
#' fmz: feature m/z
#' frt: feature RT in seconds
#' specObject: indicates if spObject contains XCMS peaks ("raw")
#' or a RAMClustR pseudo-MS/MS spectrum ("cluster")
#' mztol: absolute tolerance for feature m/z search
#' rttol: absolute tolerance for feature RT search
#' Goncalo Graca and Yuheng Cai (Rene) (Imperial College London)
#' @export

checkIsotope <- function(fmz, frt, spec, rttol = 5, mztol = 0.01){
  iso <- 0
  if(sum(abs(spec[,"mz"] - (fmz - 1.0034)) < mztol) > 0) {
    if(spec[,"into"][which.min(abs(spec[,"mz"] - (fmz-1.0034)))] >
       spec[,"into"][which.min(abs(spec[,"mz"] - fmz))]){
        iso <- 1
      if(sum(abs(spec[,"mz"] - (fmz - 1.0034 * 2)) < mztol) > 0){
        if(spec[,"into"][which.min(abs(spec[,"mz"] - (fmz-1.0034)))] >
           spec[,"into"][which.min(abs(spec[,"mz"] - (fmz - 1.0034)))]){
            iso <- 2
          if(sum(abs(spec[,"mz"] - (fmz - 1.0034 * 3)) < mztol) > 0){
            if(spec[,"into"][which.min(abs(spec[,"mz"] - (fmz - 1.0034)))] >
               spec[,"into"][which.min(abs(spec[,"mz"] - (fmz - 1.0034 * 2)))]){
              iso <- 3
            }
          }
        }
      }
    }
  }
  return(iso)
}
