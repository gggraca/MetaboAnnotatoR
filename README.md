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

## Vignettes
An example of usage using one LC-MS AIF chromatogram is provided in the [introductory vignette](/vignettes/introduction.html).
An illustration of the generation of Metabolite database records for MetaboAnnotatoR is given [here](/vignettes/gen_library_entry.html).
