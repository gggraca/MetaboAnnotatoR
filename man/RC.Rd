\name{RC}
\alias{RC}
\docType{data}
\title{Pseudo-MS/MS spectra from MESA Lipid Positive LC-MS serum samples}
\description{
Pseudo-MS/MS spectra obtained by processing the MESA Lipid Positive LC-MS xcms data using  
(\href{https://pubs.acs.org/doi/abs/10.1021/ac501530d}{RAMClustR package}). 
The xcms data consisted random selection of 100 LC-MS Lipid Positive chromatograms from human serum samples acquired for 
the Multi-Ethnic Study of Atherosclerosis (\href{https://mesa-nhlbi.org/}{MESA}).
}

\details{
One hundred serum samples from the MESA study, were randomly selected and the MS1 and MSE scans from Lipid Positive LC-MS
chromatograms were separated in 2 different files for processing in xcms and RAMClustR.
The raw data was acquired on a Waters Acquity UPLC system connected to a Waters Xevo-G2 Q-ToF system 
operated in the MSE mode (all-ion fragmentation mode), and processed in xcms. 
The full experimental details, xcms and RAMClustR parameter sets used for processing are described in the 
\href{https://pubs.acs.org/doi/10.1021/acs.analchem.1c03032}{MetaboAnnotatoR paper}.
}

\format{A RAMClustR object (list).}
\usage{data("RC")}
\references{
Graca G., Cai Y., Lau C-H. E., Vorkas P.A., Lewis M.R., Want E.J., Herrington D., Ebbels, T.M.D.
Automated Annotation of Untargeted All-Ion Fragmentation LC-MS Metabolomics Data with MetaboAnnotatoR. \emph{Analytical Chemistry}, 
2022, 94(8), 3446-3455. DOI: 10.1021/acs.analchem.1c03032
}
