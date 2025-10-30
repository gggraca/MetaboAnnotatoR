## test for RCspec function

# read test data
RCpath <- system.file("/Data/MESA_RAMClustR.rds", package = "MetaboAnnotatoR")
RC <- readRDS(RCpath)
expected <- data.frame(
    mz=c(391.2248, 469.3115, 468.3089, 450.2982, 471.3170, 451.3013),
    into=c(1015.122, 303973.046, 1279964.631, 36714.991, 7463.412, 8979.110),
    rt=c(81.577, 81.612, 81.827, 82.351, 82.380, 82.729))

test_that("Correct Pseudo-MS/MS spectrum is retrieved from RAMClustR object", {
  pseudoSpec <- RCspec(fmz=468.3094, frt=82.92, RC) 
  expect_equal(pseudoSpec, expected)
})