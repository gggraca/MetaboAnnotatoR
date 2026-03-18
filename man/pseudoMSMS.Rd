\name{pseudoMSMS}
\alias{pseudoMSMS}
\docType{data}
\title{MetaboAnnotatoR-generated example pseudoMSMS spectrum and related data}
\description{
This dataset is an example of a MetaboAnnotatoR-generated pseudoMSMS spectrum for an
LC-MS feature of interest (468.3095 m/z, 82.92 s) from a serum Lipidomics LC-MS dataset (Quality control sample).
The raw data was acquired on a Waters Acquity UPLC system connected to a Waters Xevo-G2 Q-ToF system 
operated in the MSE mode (all-ion fragmentation mode).
}
\details{
The raw mzML file from which the data was extracted can be downloaded from zenodo:
\href{https://zenodo.org/records/17408169/files/Lipid_Positive_QC.mzML}{Lipid_Positive_QC.mzML}.
}
\usage{data("pseudoMSMS")}
\format{
A list of 7 elements related to the feature of interest: 
1) in-source fragment pseudo-MS/MS spectrum from MS1 scans (data frame);
2) all-ion fragmentation (AIF) pseudo-MS/MS spectrum (data frame);
3) MS1 spectrum (data frame);
4) AIF spectrum (data frame);
5) AIF extracted-ion chromatograms (EICs) (MChromatograms object);
6) m/z from the AIF spectrum (vector);
7) EIC from the feature of interest (MChromatograms object).
}
\seealso{
\code{\link[MetaboAnnotatoR]{getPseudoMSMS}}
\code{\link[MSnbase]{MChromatograms}}
}
\author{
Goncalo Graca
}
