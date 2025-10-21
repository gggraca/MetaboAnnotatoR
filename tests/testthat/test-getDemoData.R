test_that("demo files successfuly moved to working directory", {
  getDemoData()
  expect_true(file.exists("targetTable.csv"))
  expect_true(file.exists("XCMS_options.csv"))
  # remove files created by function
  on.exit( tryCatch({ file.remove(c("targetTable.csv","XCMS_options.csv")) }, error=function(e){ invisible() }, warning=function(w){ invisible() }) )
})