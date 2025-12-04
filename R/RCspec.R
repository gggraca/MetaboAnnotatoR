#' Get RAMClustR pseudo-MS/MS spectra (cluster).
#'
#' For a given feature, find out corresponding cluster (pseudo-MS/MS spectra) 
#' from RAMClustR object.
#'
#' @author Goncalo Graca & Yuheng (Rene) Cai (Imperial College London)
#'
#' @param fmz The m/z for the feature of interest.
#' @param frt Retention time in seconds for the feature of interest.
#' @param ramclustObj RAMClustR object (list) with parent-fragment 
#' reconstructions. See RAMClustR paper for more details 
#' (https://pubs.acs.org/doi/10.1021/ac501530d).
#' @return Pseudo-MS/MS spectrum for the feature of interest.
#' @examples
#' # read in ramclustR object
#' RCpath <- system.file("/Data/MESA_RAMClustR.rds", package="MetaboAnnotatoR")
#' ramclustObj <- readRDS(RCpath)
#' # obtain pseudo-MS/MS
#' spec <- RCspec(fmz=468.3094, frt=82.92, ramclustObj)
#' @export
RCspec <- function(fmz, frt, ramclustObj){
    # find a given feature in RAMClustR object by mass
    pseudoSpec <- NULL
    altclus <- findFeature(ramclustObj, mz=fmz, mztol=0.02, rttol=2)
    if(length(altclus$rt) > 0){
        # select matched cluster with nearest rt
        if(min(abs(altclus$rt - frt)) < 3){
            selectclus <- which.min(abs(altclus$rt - frt))
            nclus <- altclus$featclus[selectclus] # id of matched cluster
            # check if the matched cluster is cluster 0
            if(nclus > 0){
                # extract pseudo spectra
                inclus <- which(ramclustObj$featclus == nclus)
                pseudoSpec <- data.frame("mz"=ramclustObj$fmz[inclus],
                                            "into"=ramclustObj$msint[inclus],
                                            "rt"=ramclustObj$frt[inclus])
            }
        }
    }
    return(pseudoSpec)
}

## helper function-------------------------------------------------------------
# findFeature: finds the correct feature in the RAMClustR object
findFeature <- function(ramclustObj, mz, mztol=0.02, rttol=2){
    target <- which(abs(ramclustObj$fmz-mz) <= mztol)
    if(length(target) == 0) {
        out <- data.frame(featn=NA, featclus=NA, mz=NA, rt=NA, 
                        int=NA, M0_plausible=NA)
        out<-out[0,]
    } else {
        out <- data.frame(featn=target, featclus=ramclustObj$featclus[target], 
                        mz=ramclustObj$fmz[target], rt=ramclustObj$frt[target], 
                        int=ramclustObj$msint[target], 
                        M0_plausible=rep(NA, length(tar)))
    }
    return(out)
}
