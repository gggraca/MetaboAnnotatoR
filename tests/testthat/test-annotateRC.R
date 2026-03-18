## test for annotateRC function

## load test data
data("RC")
data("xset")

## load library
data("LipidPos")

targets <- data.frame(feature.mz=520.3408533, feature.rt=100.6238759)

expected <- data.frame(feature.mz=520.3408533,
					   feature.rt=100.6238759,
					   metabolite="LPC(18:2)",
					   feature.type="parent",
					   ion.type="[M+H]+",
					   isotope="M+0",
					   mz.metabolite=520.339805,
					   matched.mz=520.339805,
					   mz.error=2.01464096744439,
					   pseudoMSMS="FALSE",
					   fraction="3  of  4",
					   score=0.423183179077441
					   )

test_that("output files are created", {
    # run function  
    suppressMessages(
    	annotations <- annotateRC(targets, 
                                  xcmsObject=xset, ramclustObj=RC, 
    							  libs="LipidPos", RTs="none", checkIsotope=TRUE)
      )
	global <- annotations$global
    # test if annotations object is a list 
    expect_true(is.list(annotations))
    
    # test if expected is equal to annotations$global result object
    expect_equal(global, expected)
})
