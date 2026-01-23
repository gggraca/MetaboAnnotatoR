## test for annotateRC function

## load test data
data("MESA_RAMClustR")
data("MESA_Xset")
tPath <- system.file("targetTable.csv", package = "MetaboAnnotatoR")

test_that("output files are created", {
    # run function  
    suppressMessages(
      annotateRC(targetTable=tPath, 
                 xcmsObject=xset, 
                 ramclustObj=RC, 
                 libs="Lipids", 
                 ESImode="POS")
      )
    # test presence of Annotation folder 
    expect_true(dir.exists("./Annotations/"))
    # test the size of the result file list
    resultsPath <- paste("./Annotations/", dir("./Annotations/"), sep="")
    fileList <- dir(resultsPath)
    expect_equal(length(fileList), 22)
    on.exit( tryCatch({ file.remove(fileList) }, 
                      error=function(e){ invisible() }, 
                      warning=function(w){ invisible() }))
})

# remove ./Libraries and files
unlink("Annotations", recursive=TRUE)
