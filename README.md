# MetaboAnnotatoR

<h3> Description </h3>
This R package is designed to perform metabolite annotation of features from LC-MS All-ion fragmentation (AIF) datasets, using ion fragment databases.
It requires raw LC-MS AIF chromatograms acquired/transformed in centroid mode or processed data outputs obtained using [RAMClustR](https://github.com/cbroeckl/RAMClustR).

<h3> Installation instructions </h3>

To install via Github, run the following code in R console:
```
install.packages("devtools", repos = "http://cran.us.r-project.org", dependencies = TRUE)

library(devtools)

install_github("gggraca/MetaboAnnotatoR", dependencies = TRUE)
```

<h3> Example session </h3>
An example of usage using one LC-MS AIF chromatogram is provided below.

First load the package and dependencies:
```
library(MetaboAnnotatoR)
```
As an input the software requires three files, that should be present in the working directory:

1. A .csv file containing the features to be annotated (targetTable.csv - the file name can be changed);
2. A .csv file containing the XCMS peak-picking options (XCMS_options.csv);
3. A a raw chromatogram file in .mzML format containing low and high (AIF) collision energy scans or two .netCDF files, 
each one containing the low or high collision energy scans. 

The following function will place three example files as described above in the user working directory:
```
getDemoData()
```
The dataset used consists of a lipidomics (reverse-phase chromatography) dataset acquired in ESI-MS positive mode. It contains both lipids and lipid-like molecules.

The 'targetTable.csv' file contains 6 features, corresponding to 4 phospholipids and 2 lipid-like molecules.

Let's start the annotation by searching the Lipid libraries:

```
annotateAIF(targetTable = "targetTable.csv", 
  filetype = "mzML", 
  libs = "Lipids",
  ESImode = "POS",
  RTfile = "none",
  nCE = 1,
  corThresh = 0.7,
  checkIsotope = TRUE)
```
A new folder will be automatically created in the working directory called  'Annotation', were annotation results will be stored.
This folder contains several files. Annotation results are stored in .csv files for each feature. The annotated 'targetTable' 
can be found in the 'Global_Results.csv' file. Graphical results are stored as .pdf files in the same folder.

Some fetures will remain unannotated, which means that no matching compound was found in the lipids library. 

We can now try to use the small molecule libraries to annotate :

```
annotateAIF(targetTable = "targetTable.csv", 
  filetype = "mzML", 
  libs = "Metabolites",
  ESImode = "POS",
  RTfile = "none",
  nCE = 1,
  corThresh = 0.7,
  checkIsotope = TRUE)
```
The libraries used in the example are the default 'Lipids' and 'Metabolites' libraries provided with the package.
