## test for annotateRC function

## input test data
fpath <- system.file("/Data/MESA_RAMClustR.rds", package = "MetaboAnnotatoR")
RC <- readRDS(fpath)

# create new xcmsSet and add a matrix with peaks
# note that ms1 and aif (ms2) peaks must be added so that xcmsSpec function works properly
# "dummy" fragments were added to the Xset to avoid having in Source Spectra with one ion only 
Xset <- new("xcmsSet")
tPath <- system.file("targetTable.csv", package = "MetaboAnnotatoR")
ms1 <- read.csv(tPath)
colnames(ms1) <- c("mz", "rt", "into")
ms1$into <- 1
ms1$sample <- 1
dummy <- ms1
dummy$into <- dummy$into/2
dummy$mz <- dummy$mz/2
ms2 <- ms1
ms2$sample <- 2
dummyFrags <- ms2
dummyFrags$mz <- dummyFrags$mz/2
dummyFrags$into <- dummyFrags$into/2
Xset@peaks <- as.matrix(rbind(ms1, ms2, dummy, dummyFrags))

test_that("output files are created", {
    # run function  
    suppressMessages(
      annotateRC(targetTable=tPath, xcmsObject=Xset, ramclustObj=RC, libs="Lipids", ESImode="POS")
      )
    # test presence of Annotation folder 
    expect_true(dir.exists("./Annotations/"))
    # test the size of the result file list
    resultsPath <- paste("./Annotations/", dir("./Annotations/"), sep="")
    fileList <- dir(resultsPath)
    expect_equal(length(fileList), 18)
    on.exit( tryCatch({ file.remove(fileList) }, 
                      error=function(e){ invisible() }, warning=function(w){ invisible() }))
})

# remove ./Libraries and files
unlink("Annotations", recursive=TRUE)
