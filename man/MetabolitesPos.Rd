\name{MetabolitesPos}
\alias{MetabolitesPos}
\docType{data}
\title{MetaboAnnotatoR Metabolite Positive fragment libraries}
\description{
MetaboAnnotatoR default nonlipid small-molecule positive electrospray ionisation mode (ESI+) fragment libraries, 
comprising 102 metabolite class/adducts combinations corresponding to a total of 255 entries. 
The libraries consist of records that include parent ion m/z and expected fragments from collision-induced decay (CID) MS/MS experiments. 
The nonlipid small-molecule ESI+ library was generated from experimental CID MS/MS spectra from proton or sodium adduct ions 
corresponding to metabolites commonly found in human biofluids, such as urine and blood serum or plasma deposited in 
\href{https://massbank.eu/}{MassBank} and \href{https://gnps.ucsd.edu/ProteoSAFe/libraries.jsp}{GNPS} databases.
}
\details{
Each fragment entry from contains the fragment m/z values of each metabolite and respective score weights, 
which were attributed based on the experimental intensity of the reference MS/MS spectra from MassBank and GNPS. 
The full details on the libraries are described in the  
\href{https://pubs.acs.org/doi/10.1021/acs.analchem.1c03032}{MetaboAnnotatoR paper}.
This object contains both the libraries for the different metabolite classes/adducts and the filepaths of the original
libraries .csv files.
}
\usage{MetabolitesPos}
\format{A list containing two large lists: 1) one library (data frame) per metabolite/adduct combination; 
2) A vector with the original libraries filenames.}

\references{
Graca G., Cai Y., Lau C-H. E., Vorkas P.A., Lewis M.R., Want E.J., Herrington D., Ebbels, T.M.D.
Automated Annotation of Untargeted All-Ion Fragmentation LC-MS Metabolomics Data with MetaboAnnotatoR. \emph{Analytical Chemistry}, 
2022, 94(8), 3446-3455. DOI: 10.1021/acs.analchem.1c03032
}
