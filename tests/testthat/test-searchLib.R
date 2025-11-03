# expected result data
fmz <- 152.0720
frt <- 125
filePath <- system.file("/Data/Acetaminophen_pos.csv", package = "MetaboAnnotatoR")
expected <- list()
expected[[1]] <- read.csv(filePath, header = TRUE, sep=",", check.names = FALSE)

# test data
filePath2 <- system.file("/Data/Metabolites_POS.rds", package = "MetaboAnnotatoR")
libraries <- readRDS(filePath2)
libfiles <- libraries$libfiles
lib <- libraries$lib
inSourceSpec <- data.frame(mz = 152.0720, into = 1)

test_that("Obtains the correct candidates", {
  result <- searchLib(lib, libfiles, fmz, frt, tolerance = 25, RTs = "none", inSourceSpec)
  expect_true(is.list(result))
  expect_equal(result[[1]], expected[[1]])
})
