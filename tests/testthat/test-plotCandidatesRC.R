# input test data
fmz <- 152.0723
frt <- 125
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

highCESpec <- data.frame(mz=c(59.0489, 65.0389, 66.0427, 67.0550, 70.0659, 73.0762, 82.0658, 92.0498, 93.0355, 93.0569, 
                              109.0523, 110.0622, 111.0452, 111.0647, 112.0476, 121.0408, 134.0611, 136.0762, 152.0716, 154.0781), 
                         into=c(3228, 8696, 564, 1004, 432, 592, 2092, 4836, 832, 560, 448, 30696, 8516, 3400, 464, 804, 4480, 368, 65236, 464))

# run rankScore function and get the ranked result
rankedCandidates <- rankScore(result, specMatch)

test_that("plot successfuly saved to file", {
  plotCandidatesRC(fmz, frt, highCESpec, DatasetName = "test_dataset",
                               rankedCandidates, candidate=1, DirPath="./")
  expect_true(file.exists("test_dataset_152.072mz_125_candidate_1.pdf"))
  # remove files created by function
  on.exit( tryCatch({ file.remove(c("test_dataset_152.072mz_125_candidate_1.pdf")) }, 
                    error=function(e){ invisible() }, warning=function(w){ invisible() }))
})
