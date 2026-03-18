## test mspToLib

# input test data
# use tempdir as test target folder
msp_path <- system.file("extdata", "MassBank_example.msp", package = "MetaboAnnotatoR")

testLib <- tempdir()

mspToLib(msp_path,
         LibDir = testLib,
         noise = 0.005,
         mpeaksScore = 0.9,
         mpeaksThres = 0.1)

# list the files created under testLib
files <- list.files(testLib, pattern=".csv")

# some file names to check
someExpectedFiles <- c("Aflatoxin B1_positive.csv", "Aflatoxin G2_positive.csv", "Aflatoxin M1_positive.csv", "Alternaric acid_positive.csv")

test_that("msp records successfuly imported", {
  expect_contains(files, someExpectedFiles)
})