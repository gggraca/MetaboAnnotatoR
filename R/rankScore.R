#' @title Rank candidate metabolite annotations by score
#'
#' @description
#' Ranks the annotation results by final match score.
#'
#' @author Goncalo Graca & Yuheng (Rene) Cai (Imperial College London)
#'
#' @param result Results from fragment matching as data frame.
#' @param specMatch Pseudo-MS/MS of ions matched to candidate fragments.
#' @return Ranked results annotation table and ranked matched spectra as list.
#' @examples
#' # Create a result object for 3 hypothetical metabolite matches
#' result <- data.frame(metabolite=c("Met A", "Met B", "Met C"), 
#'                     feature.type=rep("parent",3), 
#'                     ion.type=rep("[M+H]+"),
#'                     isotope=rep("M+0",3),
#'                     mz.metabolite=rep(152.0723, 3), 
#'                     matched.mz=rep(152.0706, 3),
#'                     mz.error=rep(11, 3),
#'                     pseudoMSMS=rep("TRUE", 3),
#'                     fraction=c("2 of 5", "4 of 5","3 of 5"), 
#'                     score=c(0.4, 0.9, 0.6)) 
#' # Create a list of matched spectra for the 3 hypothetical metabolite matches
#' specMatch <- list()
#' specMatch$`Met A` <- data.frame(
#'                      mz=c(152.0716, 134.0611, 59.0489, 65.0389, 66.0427),
#'                      into=c(432, 592, 2092, 4836, 832)
#'                      )
#' specMatch$`Met B` <- data.frame(
#'                      mz=c(152.0716, 134.0611, 110.0622, 109.0523, 59.0489), 
#'                      into=c(65236, 4480, 30696, 448, 432)
#'                      ) 
#' specMatch$`Met C` <- data.frame(
#'                      mz=c(152.0716, 134.0611, 110.0622, 109.0523, 93.0569), 
#'                      into=c(65236, 4480, 30696, 464, 804)
#'                      )
#' # Rank candidate metabolite result and spectraMatch by score 
#' ranked <- rankScore(result, specMatch)
#' @export
rankScore <- function(result, specMatch){
    if(!is.null(result)){
        if(nrow(result)>1) {
            rankedResult <- result
            rankedResult["rank"] <- NA
        # rank
        orderScores <- order(rankedResult$score,decreasing = TRUE)
        rankedResult <-  rankedResult[orderScores,]
        rankedSpec <- specMatch[orderScores] # order spectra
        rankedResult$rank <- c(1,rep(NA,dim(rankedResult)[1]-1))
        for(i in 2:dim(rankedResult)[1]){
            if(rankedResult$score[i]==rankedResult$score[i-1]) {
                rankedResult$rank[i]<-rankedResult$rank[i-1]
            } else if(rankedResult$score[i]<rankedResult$score[i-1]) {
                rankedResult$rank[i]<-rankedResult$rank[i-1]+1
            }
        }
        } else if (nrow(result)==1) {
            rankedResult <- result
            rankedResult["rank"] <- 1
            rankedSpec <- specMatch
        }
    }
    if (is.null(result)){
        rankedResult <- NULL
        rankedSpec <- NULL
        }
    return(list("rankedResult"=rankedResult,"rankedSpecMatch"=rankedSpec))
}
