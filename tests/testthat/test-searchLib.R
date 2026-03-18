## test searchLib function

# expected result data
fmz <- 152.0720
frt <- 125

# load library data
data("MetabolitesPos")

# load acetaminophen data
data("acetaminophen")

lib <- MetabolitesPos$lib
libfiles <- MetabolitesPos$libfiles
inSourceSpec <- data.frame(mz = 152.0720, into = 1)

test_that("Obtains the correct candidates", {
    result <- searchLib(lib, libfiles, fmz, frt, tolerance = 25, RTs = "none", inSourceSpec)
    expect_true(is.list(result))
    expect_equal(result[[1]], acetaminophen[[1]])
})
