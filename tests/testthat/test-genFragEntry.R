# read expected table
testFile <- system.file("/Data/Metabolites_POS.rds", package = "MetaboAnnotatoR")
expected <- readRDS(testFile)
expected <- expected$lib[[69]]

# input test spectrum as data frame
specObject <- data.frame(V1=c(70.0298, 85.0652, 90.0556, 98.024, 116.0353, 124.0766, 184.0981, 202.1085, 220.1185),
                         V2=c(13.965907, 13.534607, 100.0, 26.165537, 15.383036, 25.231054, 28.578764, 43.017047, 64.962005))

test_that("library entry correctly generated", {

  entry <- genFragEntry(specObject, 
                        "Pantothenic acid",
                        "[M+H]+",220.1179,
                        "Pantothenic_acid_pos", 
                        noise=0.005, 
                        mpeaksScore=0.9, 
                        mpeaksThres=0.1, 
                        mzTol=0.01)
  
  expect_true(file.exists("Pantothenic_acid_pos.csv"))
  
  entry <- read.csv("Pantothenic_acid_pos.csv", check.names=FALSE)
  
  expect_equal(entry, expected)
  
  # remove file created by function
  on.exit( tryCatch({ file.remove('./Pantothenic_acid_pos.csv') }, error=function(e){ invisible() }, warning=function(w){ invisible() }) )
})
