# read test data
spObject <- data.frame(
  mz = c(703.5769, 704.5799, 705.5812, 706.5721),
  into = c(205458624, 85536216, 22717336, 5887723)
  )

test_that("correct isotope label is detected", {
  expected_iso <- 0
  obtained_iso <- checkIsotope(fmz = 703.5769, 
                               frt = 70, 
                               spec = spObject,
                               rttol = 5, 
                               mztol = 0.01)
  expect_equal(expected_iso, obtained_iso)
})