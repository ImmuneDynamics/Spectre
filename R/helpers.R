#' Check packages are installed
#' 
#' To check required packages are installed
#' 
check_packages_installed <- function(packages_name) {
    for (p in packages_name) {
        if (!is.element(p, installed.packages()[, 1])) 
            stop(paste(p, "is required but not installed"))
    }
}
