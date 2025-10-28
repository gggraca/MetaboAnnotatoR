# expected result data
fmz <- 152.0720
frt <- 125
file <- system.file("/Libraries/Metabolites/POS/Acetaminophen_pos.csv", package = "MetaboAnnotatoR")
expected <- list()
expected[[1]] <- read.csv(file, header = TRUE, sep=",", check.names = FALSE)

# test data
LibPath <- system.file("/Libraries/Metabolites/POS/", package = "MetaboAnnotatoR")
libfiles <- list.files(LibPath, full.names = TRUE)
lib <- lapply(libfiles, read.csv, header = TRUE, sep=",", check.names = FALSE)
inSourceSpec <- data.frame(mz = 152.0720, into = 1)

test_that("obtains the correct candidates", {
  result <- searchLib(lib, libfiles, fmz, frt, tolerance = 25, RTs = "none", inSourceSpec)
  expect_true(is.list(result))
  expect_equal(result[[1]], expected[[1]])
})
