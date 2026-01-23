#' @import utils
.onLoad <- function(libname = find.package("MetaboAnnotatoR"),
					pkgname = "MetaboAnnotatoR"){
	# R CMD check workaround for ggplot2 variable names
	if(getRversion() >= "2.15.1") {
		utils::globalVariables(c("mz","rt","into"))
		invisible()
	}
}