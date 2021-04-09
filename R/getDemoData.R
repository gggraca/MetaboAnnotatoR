#' Function to place demonstration data into working directory.
#'
#' Places example raw chromatogram (mzML), XCMS_options and target features
#' files into the working directory to test automated AIF annotation using
#' annotateAIF function.
#'
#' @author Goncalo Graca (Imperial College London)
#'
#' @return Files: Lipid_Positive_QC.mzML, targetTable.csv, XCMS_options.csv
#' @export
getDemoData <- function() {
  chromPath <- system.file("/Data/Lipid_Positive_QC.mzML",
                           package = "MetaboAnnotatoR")
  targetTablePath <- system.file("targetTable.csv",
                                 package = "MetaboAnnotatoR")
  xcmsOptionsPath <- system.file("XCMS_options.csv",
                                 package = "MetaboAnnotatoR")

  file.copy(from = chromPath, to = getwd())
  file.copy(from = xcmsOptionsPath, to = getwd())
  file.copy(from = targetTablePath, to = getwd())
}
