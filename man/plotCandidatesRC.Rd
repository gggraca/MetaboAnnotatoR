% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotCandidatesRC.R
\name{plotCandidatesRC}
\alias{plotCandidatesRC}
\title{Visual results of annotations based on LC-MS AIF processed using RAMClusteR.}
\usage{
plotCandidatesRC(
  fmz,
  frt,
  highCESpec,
  DatasetName,
  rankedCandidates,
  candidate,
  DirPath
)
}
\arguments{
\item{fmz}{The m/z for the feature of interest.}

\item{frt}{Retention time in seconds for the feature of interest.}

\item{highCESpec}{MS2 peaks at the RT window of the feature of interest.}

\item{DatasetName}{Name of the data set processed using RAMClustR.}

\item{rankedCandidates}{List containing the ranked candidate annotations.}

\item{candidate}{Number of candidate annotations to plot.}

\item{DirPath}{Path to the folder where the plots will be saved.}
}
\value{
Saves the plots of the candidate annotations as one pdf file:
a plot pseudo-MS/MS spectrum for the matched ions.
}
\description{
Plots pseudo-MS/MS spectra of matched fragments from RC objects.
Runs inside annotateRC function.
}
\author{
Goncalo Graca & Yuheng (Rene) Cai (Imperial College London)
}
