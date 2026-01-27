## test for xcmsSpec function

# feature to test
fmz <- 585.2692
frt <- 72.8

# input expected result from xset data

expected <- xset@peaks[which(xset@peaks[,"sample"] == 7),]
idx <- which(expected[,"rt"] > (frt - 1) &
               expected[,"rt"] < (frt + 1))
expected <- expected[idx, c("mz","into")]

test_that("xcmsSpec is correctly imported", {
    result <- xcmsSpec(fmz, frt, mztol=0.01, xset, rttol=1, highCE = FALSE)
    expect_equal(result, expected)
})
