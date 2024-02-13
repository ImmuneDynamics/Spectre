#' Check packages are installed
#' 
#' To check required packages are installed
#' 
check_packages_installed <- function(packages_name) {
    for (p in packages_name) {
        if (!requireNamespace(p, quietly = TRUE)) {
            # call. set to FALSE so we don't show the code that test the package's 
            # availability
            stop(paste(p, "is required but not installed"), call.=FALSE)
        }
    }
}
