# input test data
# folders containing imported files will be created at:
# ./Libraries/Custom/POS and ./Libraries/Custom/NEG
msp_path <- system.file("/Data/MassBank_example.msp", package = "MetaboAnnotatoR")
mspToLib(msp_path,
         library_name = "Custom",
         noise = 0.005,
         mpeaksScore = 0.9,
         mpeaksThres = 0.1)

# list the files created under ./Libraries/Custom/POS
files <- list.files("./Libraries/Custom/POS")

# some file names to check
someExpectedFiles <- c("Aflatoxin B1.csv","Isomarticin.csv","Penicillic Acid.csv")

test_that("msp records successfuly imported", {
  expect_contains(files, someExpectedFiles)
})

# remove ./Libraries and files
unlink("Libraries", recursive=TRUE)