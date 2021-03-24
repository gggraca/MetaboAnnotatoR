## Function RCspec() ----------
## For a given feature, find out corresponding cluster from RAMClust and extract pseudo spectra
### Goncalo Graca & Yuheng Cai (Imperial College London)
### g.gomes-da-graca@imperial.ac.uk

RCspec <- function(fmz, frt, ramclustObj){
  library(RAMClustR)
  # find a given feature in RAMClust by mass
  pseudoSpec <- NULL
  altclus <- findmass(ramclustObj, mz = fmz, mztol = 0.02, rttol = 2, zmax = 6, m.check = TRUE)
  if (length(altclus$rt) > 0){
    # select matched cluster with nearest rt
    if (min(abs(altclus$rt - frt)) < 3){
      selectclus <- which.min(abs(altclus$rt - frt))
      nclus <- altclus$featclus[selectclus] # id of matched cluster
      # check if the matched cluster is cluster 0
      if (nclus > 0){
        # extract pseudo spectra
        inclus <- which(RC$featclus == nclus)
        pseudoSpec <- data.frame("mz" = RC$fmz[inclus], "into" = RC$msint[inclus], "rt" = RC$frt[inclus])
      }
    }
  }
  return(pseudoSpec)
}