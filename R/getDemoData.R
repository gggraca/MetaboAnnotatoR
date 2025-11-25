#' Function to place demonstration data into working directory.
#'
#' Places XCMS_options and target features files into the working directory 
#' to test automated AIF annotation using annotateAIF function.
#'
#' @author Goncalo Graca (Imperial College London)
#' @param AIFdata If AIF test data is to be loaded. Default is FALSE.
#' @return Files: targetTable.csv, XCMS_options.csv and two additional 
#' environment objects if AIFdata=TRUE: Xset (xcmsSet object) and 
#' RC (RAMClustR object).
#' @examples
#' getDemoData()
#' @export
getDemoData <- function(AIFdata=FALSE) {
    targetTablePath <- system.file("targetTable.csv", 
                                    package="MetaboAnnotatoR")
    xcmsOptionsPath <- system.file("XCMS_options.csv", 
                                    package="MetaboAnnotatoR")
    file.copy(from=xcmsOptionsPath, to=getwd())
    file.copy(from=targetTablePath, to=getwd())
    if(isTRUE(AIFdata)){
        rcpath <- system.file("/Data/MESA_RAMClustR.rds", 
                                package="MetaboAnnotatoR")
        RC <- readRDS(rcpath)
        XsetPath <- system.file("/Data/MESA_Xset.rds", 
                                package="MetaboAnnotatoR")
        Xset <- readRDS(XsetPath)
        data <- list(RC=RC, Xset=Xset)
        list2env(data, envir=.GlobalEnv)
    }
}
