#' Isotopologue type check.
#'
#' Checks the type of isotope of the input feature.
#'
#' @author Goncalo Graca and Yuheng (Rene) Cai (Imperial College London)
#'
#' @param fmz Feature m/z.
#' @param frt Feature RT in seconds.
#' @param specObject Indicates if spObject contains XCMS peaks ("raw") or a
#' RAMClustR pseudo-MS/MS spectrum ("cluster").
#' @param mztol Absolute tolerance for feature m/z search in Da.
#' @param rttol Absolute tolerance for feature RT search in seconds.
#' @return A "tag" of the isotope from an isotopic series as 0, 1, 2 or 3 for
#' M+0, M+1, M+2 and M+3, respectively.
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
