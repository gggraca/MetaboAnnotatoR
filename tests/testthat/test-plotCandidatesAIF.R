## test for plotCandidatesAIF function

## input test data
fmz <- 468.3094542
frt <- 82.92008607

highCESpec <- pseudoMSMS$ms2
ms2eic <- pseudoMSMS$ms2_eic
rankedCandidates <- rCandidates

test_that("plot successfuly saved to file", {
    plotCandidatesAIF(fmz, frt, highCESpec, ms2eic, SpName="AIF_test", rankedCandidates, candidate=1, DirPath="./")
    expect_true(file.exists("AIF_test_468.309mz_82.92_candidate_1.pdf"))
    ## remove files created by function
    on.exit( tryCatch({ file.remove(c("AIF_test_468.309mz_82.92_candidate_1.pdf")) }, 
                      error=function(e){ invisible() }, warning=function(w){ invisible() }))
})
