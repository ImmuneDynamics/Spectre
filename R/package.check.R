#' package.check - a function to check the installation of all required packages.
#'
#' This function allows you to check to see if all the common use packages dependencies for Spectre are installed.
#' See \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
#'
#' @return returns an error message if one of the common use packages are not installed. Proceeds in order of package importance, and only the first error message encountered will be returned.
#'
#' @param type DEFAULT = "general". If "general", then checks for the packages required for general Spectre usage. If "spatial", then checks for additional packages required for spatial analysis. If "ML", then checks for additional packages required for machine-learing functionality.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' 
#'
#' @examples
#' package.check()
#'
#' @export

package.check <- function(type = "general") {
  
    deets <- utils::packageDescription("Spectre")
    os.deets <- utils::sessionInfo()
    r.deets <- R.Version()
    
    message(paste0("Package: ", deets$Package))
    message(paste0(" -- Version:           ", deets$Version))
    message(paste0(" -- Install date:      ", deets$Packaged))
    message(paste0(" -- Install source:    ", deets$RemoteType))
    message(paste0(" -- R version:         ", r.deets$version.string))
    message(paste0(" -- OS:                ", os.deets$running))
    message(paste0(" -- OS detail:         ", os.deets$platform))
    message(" -- Library path(s):      ")
    
    for(i in .libPaths()){
      message(paste0("        ", i))
    }

    message(paste0("               "))
    
    # if(deets$ondiskversion != deets$loadedversion){
    #   message(paste0("Please note, your "on disk" version of Spectre is different to the "loaded" version (i.e. you appear to have had the Spectre library loaded and active while installing a new version. To address this, either:"))
    #   message(paste0("    a) run "detach("package:Spectre", unload = TRUE)" followed by "library("Spectre"), or"))
    #   message(paste0("    b) re-start RStudio"))
    #   message("   ")
    # }
    
    message(paste0(("Checking dependency packages...")))

    imports_split <- strsplit(deets$Imports, ",")[[1]]
    # Remove version constraints and trim whitespace
    imports_clean <- trimws(gsub("\\s*\\(.*?\\)", "", imports_split))
    
    res.list <- lapply(imports_clean, function(x) {
        if (!is.element(x, installed.packages()[,1])){
            return(x)
        } else {
            return(NULL)
        }
    })
    res.list <- unlist(res.list)
  
    ### Final message
    
    if (length(res.list) > 0) {
        warning(paste0("Please install the following packages: ", paste(res.list, collapse = ", ")))
        message("Check out https://immunedynamics.github.io/spectre/getting-started/ for help with installation")
    } else {
        message(" -- All required packages are installed.")
    }

    if (type == "spatial"){
        cat("\n\n")
        cat("Checking additional packages required for spatial analysis...")
        spatial.packages <- c("raster", "tiff", "exactextractr", "sp", "sf", "stars", "qs", "s2", "rhdf5")
        
        res.list <- lapply(spatial.packages, function(x) {
            if (!is.element(x, installed.packages()[,1])){
                return(x)
            } else {
                return(NULL)
            }
        })
        res.list <- unlist(res.list)
        
        ### Final message
        
        if (length(res.list) > 0) {
            warning(paste0("\nPlease install the following packages for spatial analysis: ", paste(res.list, collapse = ", ")))
            message("\nCheck out https://immunedynamics.github.io/spectre/getting-started/ for help with installation")
        } else {
            message("\n -- All required packages for spatial analysis are installed.")
        }
    }

}

