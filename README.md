# MetaboAnnotatoR

## Description
This R package is designed to perform metabolite annotation of features from LC-MS All-ion fragmentation (AIF) datasets, using ion fragment databases.
It requires raw LC-MS AIF chromatograms acquired/transformed in centroid mode or processed data outputs obtained using [RAMClustR](https://github.com/cbroeckl/RAMClustR).

## Installation instructions

To install via Github, run the following code in R console:
```
install.packages("devtools", repos = "http://cran.us.r-project.org", dependencies = TRUE)

library(devtools)

install_github("gggraca/MetaboAnnotatoR", dependencies = TRUE)
```
## Installation issues

In some environments, some errors might occur durring installation due to some of the package dependencies, being the most common "mzR has been built against a different Rcpp version". This issue can be easily fixed as suggested [here](https://support.bioconductor.org/p/134630/). The following environment variable can be changed before installing using install_github:
```
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
```

## Vignettes
An example of usage using one LC-MS AIF chromatogram is provided in the [introductory vignette](http://htmlpreview.github.io/?https://github.com/gggraca/MetaboAnnotatoR/blob/master/vignettes/introduction.html).
An illustration of the generation of Metabolite database records for MetaboAnnotatoR from text files is given [here](http://htmlpreview.github.io/?https://github.com/gggraca/MetaboAnnotatoR/blob/master/vignettes/gen_library_entry.html). MS/MS spectra can also be imported from .msp files as exemplified [here](http://htmlpreview.github.io/?https://github.com/gggraca/MetaboAnnotatoR/blob/master/vignettes/import_from_msp.html).
