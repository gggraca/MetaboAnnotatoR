#' Function to place demonstration data into working directory.
#'
#' Places XCMS_options and target features files into the working directory 
#' to test automated AIF annotation using annotateAIF or annotateRC functions.
#'
#' @author Goncalo Graca (Imperial College London)
#' @return Files: targetTable.csv, XCMS_options.csv
#' @examples
#' # set working directory
#' setwd(tempdir())
#' getDemoData()
#' @export
getDemoData <- function() {
    targetTablePath <- system.file("targetTable.csv", 
                                    package="MetaboAnnotatoR")
    xcmsOptionsPath <- system.file("XCMS_options.csv", 
                                    package="MetaboAnnotatoR")
    file.copy(from=xcmsOptionsPath, to=getwd())
    file.copy(from=targetTablePath, to=getwd())
}
