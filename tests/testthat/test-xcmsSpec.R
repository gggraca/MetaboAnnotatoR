## test for xcmsSpec function

# load test data from the msdata package
library(msdata)
data("xs")

# input expected result
expected <- xs@peaks[which(xs@peaks[,"sample"] == 6),c("mz","into")]


test_that("multiplication works", {
    result <- xcmsSpec(410.1440, -1, mztol=0.001, xs, rttol=1, highCE = FALSE)
    expect_equal(result, expected)
})
