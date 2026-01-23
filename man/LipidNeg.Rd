\name{LipidNeg}
\alias{libraries}
\docType{data}
\title{MetaboAnnotatoR Lipid Negative fragment libraries}
\description{
MetaboAnnotatoR default Lipid Negative fragment libraries, comprising 44 lipid class/adducts combinations,
corresponding to a total of 15729 lipid entries. The libraries consist of records that include parent ion m/z 
and expected fragments from negative electrospray ionisation (ESI) collision-induced decay (CID) MS/MS experiments.
The lipid libraries were adapted from those of the \href{https://link.springer.com/article/10.1186/s12859-017-1744-3}{LipidMatch R package} 
which is a library of theoretical m/z values for experimentally observed lipid fragments. 
The libraries were adapted to retain only fragments that were commonly observed experimentally in ESI MS/MS spectra and well-documented in the literature.
}
\details{
Each fragment entry from contains the m/z values of each lipid fragment and their respective score weights, 
which were attributed based on experimental intensity observed. The full details on the libraries are described in the  
\href{https://pubs.acs.org/doi/10.1021/acs.analchem.1c03032}{MetaboAnnotatoR paper}.
This object contains both the libraries for the different lipid classes/adducts and the filepaths of the original
libraries .csv files.
}

\format{A list containing two large lists: 1) one library (data frame) per lipid class; 2) A vector with the original libraries filenames.}
\usage{data(LipidNeg)}
\references{
Graca G., Cai Y., Lau C-H. E., Vorkas P.A., Lewis M.R., Want E.J., Herrington D., Ebbels, T.M.D.
Automated Annotation of Untargeted All-Ion Fragmentation LC-MS Metabolomics Data with MetaboAnnotatoR. \emph{Analytical Chemistry}, 
2022, 94(8), 3446-3455. DOI: 10.1021/acs.analchem.1c03032
}
