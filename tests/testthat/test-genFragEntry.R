context("genFragEntry()")

# remove file created by function
on.exit( tryCatch({ file.remove('./Pantothenic_acid_pos.csv') }, error=function(e){ invisible() }, warning=function(w){ invisible() }) )

testFile <- system.file("./Libraries/Metabolites/POS/Pantothenic acid_pos.csv",
                        package = "MetaboAnnotatoR")
expected <- read.csv(testFile, check.names = FALSE)

test_that("library entry correctly generated", {
  
  spec <- system.file("/Data/Pantothenic_acid_pos.txt", package = "MetaboAnnotatoR")
  specObject <- read.table(spec, header=FALSE)
  entry <- genFragEntry(specObject, 
                        "Pantothenic acid",
                        "[M+H]+",220.1179,
                        "Pantothenic_acid_pos", 
                        noise = 0.005, 
                        mpeaksScore = 0.9, 
                        mpeaksThres = 0.1, 
                        mzTol = 0.01)
  
  expect_true(file.exists("Pantothenic_acid_pos.csv"))
  
  entry <- read.csv("Pantothenic_acid_pos.csv", check.names = FALSE)
  
  expect_equal(entry, expected)
})