# load test data
fmz <- 152.0720
frt <- 125
iso <- 0
data("MetabolitesPos")
lib <- libraries$lib

pseudoSpec <- data.frame(mz=c(59.0489, 65.0389, 66.0427, 67.0550, 70.0659, 73.0762, 82.0658, 92.0498, 93.0355, 93.0569, 109.0523,
                                110.0622, 111.0452, 111.0647, 112.0476, 121.0408, 134.0611, 136.0762, 152.0716, 154.0781), 
                         into=c(3228, 8696, 564, 1004, 432, 592, 2092, 4836, 832, 560, 448, 30696, 8516, 3400, 464, 804, 4480, 368, 65236, 464))
highCESpec <- pseudoSpec

candidate <- lib[[8]]

# prepare expected result
expected <- list()
expected$result <- data.frame(feature.mz=fmz, 
                              feature.rt=frt, 
                              metabolite="Acetaminophen",
                              feature.type="parent",
                              ion.type="[M+H]+",
                              isotope="M+0",
                              mz.metabolite=152.0706,
                              matched.mz=152.0706,
                              mz.error=9.20616550,
                              pseudoMSMS="TRUE",
                              fraction="4  of  5",
                              score=0.529311429)
expected$specMatch <- list()
expected$specMatch$Acetaminophen <- data.frame(mz=c(152.0716, 134.0611, 110.0622, 109.0523),
                                               into=c(65236, 4480, 30696, 448))

# run function and get result
result <- compFrag(candidate, lib, fmz, frt, iso, highCESpec, pseudoSpec,
                   maxMZdiff = 0.01,
                   matchWeight = 0.5,
                   useMZerrorWeight = TRUE,
                   NoMatchWeight = 0.5,
                   additional = TRUE)

# test outputs
test_that("matching is correct", {
  expect_equal(result$result,expected$result)
})


test_that("matched fragments correct", {
  expect_equal(result$specMatch,expected$specMatch)
})