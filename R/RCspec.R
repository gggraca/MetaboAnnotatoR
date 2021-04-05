#' Get RAMClustR pseudo-MS/MS spectra (cluster).
#'
#' For a given feature, find out corresponding cluster from RAMClustR object
#' and extract pseudo-MS/MS spectra.
#'
#' @author Goncalo Graca & Yuheng (Rene) Cai (Imperial College London)
#'
#' @param fmz The m/z for the feature of interest.
#' @param frt Retention time in seconds for the feature of interest.
#' @param ramclustObj RAMClustR object with parent-fragment reconstructions
#' @return Pseudo-MS/MS spectrum for the feature of interest.
#' @export
RCspec <- function(fmz, frt, ramclustObj){
  # find a given feature in RAMClust by mass
  pseudoSpec <- NULL
  altclus <- RAMClustR::findmass(ramclustObj,
                                 mz = fmz,
                                 mztol = 0.02,
                                 rttol = 2,
                                 zmax = 6,
                                 m.check = TRUE)
  if (length(altclus$rt) > 0){
    # select matched cluster with nearest rt
    if (min(abs(altclus$rt - frt)) < 3){
      selectclus <- which.min(abs(altclus$rt - frt))
      nclus <- altclus$featclus[selectclus] # id of matched cluster
      # check if the matched cluster is cluster 0
      if (nclus > 0){
        # extract pseudo spectra
        inclus <- which(RC$featclus == nclus)
        pseudoSpec <- data.frame("mz" = RC$fmz[inclus],
                                 "into" = RC$msint[inclus],
                                 "rt" = RC$frt[inclus])
      }
    }
  }
  return(pseudoSpec)
}
