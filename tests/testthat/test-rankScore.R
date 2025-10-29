# input test data
result <- data.frame(metabolite=c("Met A", "Met B", "Met C"), 
                 feature.type=rep("parent",3), 
                 ion.type=rep("[M+H]+"),
                 isotope=rep("M+0",3),
                 mz.metabolite=rep(152.0723, 3), 
                 matched.mz=rep(152.0706, 3),
                 mz.error=rep(11, 3),
                 pseudoMSMS=rep("TRUE", 3),
                 fraction=c("2 of 5", "4 of 5","3 of 5"), 
                 score=c(0.4, 0.9, 0.6)) 

specMatch <- list()
specMatch$`Met A` <- data.frame(mz=c(152.0716, 134.0611, 59.0489, 65.0389, 66.0427), into=c(432, 592, 2092, 4836, 832))
specMatch$`Met B` <- data.frame(mz=c(152.0716, 134.0611, 110.0622, 109.0523, 59.0489), into=c(65236, 4480, 30696, 448, 432)) 
specMatch$`Met C` <- data.frame(mz=c(152.0716, 134.0611, 110.0622, 109.0523, 93.0569), into=c(65236, 4480, 30696, 464, 804))

# expected results
expected <- list()
expected$result <- result[order(result$score, decreasing=TRUE),]
expected$result$rank <- 1:3
expected$specMatch <- specMatch[order(result$score, decreasing=TRUE)]

# run function and get result
ranked <- rankScore(result, specMatch)

# test outputs
test_that("results correctly ranked", {
  expect_equal(ranked$rankedResult, expected$result)
})

test_that("specMatch correctly ranked", {
  expect_equal(ranked$rankedSpecMatch, expected$specMatch)
})