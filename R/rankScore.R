## Function rankScore() ----------
## Goncalo Graca & Yuheng Cai (Imperial College London)
## g.gomes-da-graca@imperial.ac.uk

rankScore <- function(result,specMatch){
  
  if (!is.null(result)){
    if (nrow(result)>1) {
    rankedResult <- result
    rankedResult["rank"] <- NA
    
    # rank
    orderScores <- order(rankedResult$score,decreasing = TRUE)
    rankedResult <-  rankedResult[orderScores,]
    rankedSpec <- specMatch[orderScores] # order spectra
    rankedResult$rank <- c(1,rep(NA,dim(rankedResult)[1]-1))
      for(i in 2:dim(rankedResult)[1]){
        if(rankedResult$score[i]==rankedResult$score[i-1]) rankedResult$rank[i]<-rankedResult$rank[i-1]
        else if(rankedResult$score[i]<rankedResult$score[i-1]) rankedResult$rank[i]<-rankedResult$rank[i-1]+1
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