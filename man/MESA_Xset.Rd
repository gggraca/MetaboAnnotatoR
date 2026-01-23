\name{MESA_Xset}
\alias{xset}
\docType{data}
\title{MESA cohort subset Lipid Positive LC-MS samples processed using XCMS}
\description{
A random selection of 100 LC-MS Lipid Positive chromatograms from human serum samples acquired for 
the Multi-Ethnic Study of Atherosclerosis (\href{https://mesa-nhlbi.org/}{MESA}).
}
\details{
One hundred serum samples from the MESA study, were randomly selected and the corresponding MS1 and MSE scans 
the from Lipid Positive LC-MS chromatograms were separated in 2 different files for processing in xcms.
The raw data was acquired on a Waters Acquity UPLC system connected to a Waters Xevo-G2 Q-ToF system 
operated in the MSE mode (all-ion fragmentation mode). 
The full experimental details and xcms parameter sets used for processing are described in the  
\href{https://pubs.acs.org/doi/10.1021/acs.analchem.1c03032}{MetaboAnnotatoR paper}.
}
\usage{data(MESA_Xset)}
\format{An xcmsSet object containing the xcms processed data.}
\references{
Graca G., Cai Y., Lau C-H. E., Vorkas P.A., Lewis M.R., Want E.J., Herrington D., Ebbels, T.M.D.
Automated Annotation of Untargeted All-Ion Fragmentation LC-MS Metabolomics Data with MetaboAnnotatoR. \emph{Analytical Chemistry}, 
2022, 94(8), 3446-3455. DOI: 10.1021/acs.analchem.1c03032
}
\seealso{
\code{\link[xcms]{xcmsSet}}
}
