#' package.check - a function to check the installation of all required packages.
#'
#' This function allows you to check to see if all the common use packages dependencies for Spectre are installed.
#'
#' @return returns an error message if one of the common use packages are not installed. Proceeds in order of package importance, and only the first error message encountered will be returned.
#'
#' @param type DEFAULT = "general". If "general", then checks for the packages required for general Spectre usage. If "spatial", then checks for additional packages required for spatial analysis. If "ML", then checks for additional packages required for machine-learing functionality.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' @seealso See \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
#'
#' @usage package.check()
#'
#' @export

package.check <- function(type = "general") {
  # uses devtools
  # TODO update me. This is not sustainable....
  
  deets <- devtools::package_info("Spectre", include_base = FALSE, dependencies = FALSE)
  os.deets <- devtools::session_info()

  message(paste0("Package: ", deets$package))

  if (grepl("-0", deets$ondiskversion)) {
    message(paste0(" -- Version (on disk): ", deets$ondiskversion, " (development version)"))
  } else {
    message(paste0(" -- Version (on disk): ", deets$ondiskversion))
  }

  if (grepl("-0", deets$loadedversion)) {
    message(paste0(" -- Version (loaded):  ", deets$loadedversion, " (development version)"))
  } else {
    message(paste0(" -- Version (loaded):  ", deets$loadedversion))
  }

  message(paste0(" -- Install date:      ", deets$date))
  message(paste0(" -- Install source:    ", deets$source))
  message(paste0(" -- Install path:      ", deets$path))
  message(paste0(" -- R version:         ", os.deets$platform$version))
  message(paste0(" -- OS:                ", os.deets$platform$os))
  message(paste0("               "))

  if (deets$ondiskversion != deets$loadedversion) {
    message(paste0("Please note, your 'on disk' version of Spectre is different to the 'loaded' version (i.e. you appear to have had the Spectre library loaded and active while installing a new version. To address this, either:"))
    message(paste0("    a) run 'detach('package:Spectre', unload = TRUE)' followed by 'library('Spectre'), or"))
    message(paste0("    b) re-start RStudio"))
    message("   ")
  }

  message(paste0(("Checking dependency packages...")))

  res.list <- list()

  if (!is.element("data.table", installed.packages()[, 1])) {
    message(" -- data.table is required but not installed")
    res.list[["data.table"]] <- "data.table"
  }

  if (!is.element("plyr", installed.packages()[, 1])) {
    message(" -- plyr is required but not installed")
    res.list[["plyr"]] <- "plyr"
  }

  if (!is.element("dplyr", installed.packages()[, 1])) {
    message(" -- dplyr is required but not installed")
    res.list[["data.table"]] <- "data.table"
  }

  if (!is.element("tidyr", installed.packages()[, 1])) {
    message(" -- tidyr is required but not installed")
    res.list[["tidyr"]] <- "tidyr"
  }

  if (!is.element("rstudioapi", installed.packages()[, 1])) {
    message(" -- rstudioapi is required but not installed")
    res.list[["rstudioapi"]] <- "rstudioapi"
  }

  if (!is.element("Rtsne", installed.packages()[, 1])) {
    message(" -- Rtsne is required but not installed")
    res.list[["Rtsne"]] <- "Rtsne"
  }

  if (!is.element("umap", installed.packages()[, 1])) {
    message(" -- umap is required but not installed")
    res.list[["umap"]] <- "umap"
  }

  if (!is.element("reticulate", installed.packages()[, 1])) {
    message(" -- reticulate is required but not installed")
    res.list[["reticulate"]] <- "reticulate"
  }

  if (!is.element("ggplot2", installed.packages()[, 1])) {
    message(" -- ggplot2 is required but not installed")
    res.list[["ggplot2"]] <- "ggplot2"
  }

  if (!is.element("ggthemes", installed.packages()[, 1])) {
    message(" -- ggthemes is required but not installed")
    res.list[["ggthemes"]] <- "ggthemes"
  }

  if (!is.element("scales", installed.packages()[, 1])) {
    message(" -- scales is required but not installed")
    res.list[["scales"]] <- "scales"
  }

  if (!is.element("colorRamps", installed.packages()[, 1])) {
    message(" -- colorRamps is required but not installed")
    res.list[["colorRamps"]] <- "colorRamps"
  }

  if (!is.element("RColorBrewer", installed.packages()[, 1])) {
    message(" -- RColorBrewer is required but not installed")
    res.list[["RColorBrewer"]] <- "RColorBrewer"
  }

  if (!is.element("gridExtra", installed.packages()[, 1])) {
    message(" -- gridExtra is required but not installed")
    res.list[["gridExtra"]] <- "gridExtra"
  }

  if (!is.element("ggpointdensity", installed.packages()[, 1])) {
    message(" -- ggpointdensity is required but not installed")
    res.list[["ggpointdensity"]] <- "ggpointdensity"
  }

  if (!is.element("pheatmap", installed.packages()[, 1])) {
    message(" -- pheatmap is required but not installed")
    res.list[["pheatmap"]] <- "pheatmap"
  }

  if (!is.element("ggpubr", installed.packages()[, 1])) {
    message(" -- ggpubr is required but not installed")
    res.list[["ggpubr"]] <- "ggpubr"
  }

  if (!is.element("caret", installed.packages()[, 1])) {
    message(" -- caret is required for machine learning analysis but is not installed")
    res.list[["caret"]] <- "caret"
  }

  if (!is.element("class", installed.packages()[, 1])) {
    message(" -- class is required for machine learning analysis but is not installed")
    res.list[["class"]] <- "class"
  }

  if (!is.element("flowCore", installed.packages()[, 1])) {
    message(" -- flowCore is required but not installed. Please install from BioConductor.")
    res.list[["flowCore"]] <- "flowCore"
  }

  if (!is.element("Biobase", installed.packages()[, 1])) {
    message(" -- Biobase is required but not installed. Please install from BioConductor.")
    res.list[["Biobase"]] <- "Biobase"
  }

  if (!is.element("flowViz", installed.packages()[, 1])) {
    message(" -- flowViz is required but not installed. Please install from BioConductor.")
    res.list[["flowViz"]] <- "flowViz"
  }

  if (!is.element("FlowSOM", installed.packages()[, 1])) {
    message(" -- FlowSOM is required but not installed. Please install from BioConductor.")
    res.list[["FlowSOM"]] <- "FlowSOM"
  }


  if (type == "spatial") {
    if (!is.element("raster", installed.packages()[, 1])) {
      message(" -- raster is required for SPATIAL analysis but is not installed")
      res.list[["raster"]] <- "raster"
    }

    if (!is.element("tiff", installed.packages()[, 1])) {
      message(" -- tiff is required for SPATIAL analysis but is not installed")
      res.list[["tiff"]] <- "tiff"
    }

    if (!is.element("rgeos", installed.packages()[, 1])) {
      message(" -- rgeos is required for SPATIAL analysis but is not installed")
      res.list[["rgeos"]] <- "rgeos"
    }

    if (!is.element("exactextractr", installed.packages()[, 1])) {
      message(" -- exactextractr is required for SPATIAL analysis but is not installed")
      res.list[["exactextractr"]] <- "exactextractr"
    }

    if (!is.element("sp", installed.packages()[, 1])) {
      message(" -- sp is required for SPATIAL analysis but is not installed")
      res.list[["sp"]] <- "sp"
    }

    if (!is.element("sf", installed.packages()[, 1])) {
      message(" -- sf is required for SPATIAL analysis but is not installed")
      res.list[["sf"]] <- "sf"
    }

    if (!is.element("stars", installed.packages()[, 1])) {
      message(" -- stars is required for SPATIAL analysis but is not installed")
      res.list[["stars"]] <- "stars"
    }

    if (!is.element("qs", installed.packages()[, 1])) {
      message(" -- qs is required for SPATIAL analysis but is not installed")
      res.list[["qs"]] <- "qs"
    }

    if (!is.element("s2", installed.packages()[, 1])) {
      message(" -- s2 is required for SPATIAL analysis but is not installed")
      res.list[["s2"]] <- "s2"
    }

    if (!is.element("rhdf5", installed.packages()[, 1])) {
      message(" -- rhdf5 is required for SPATIAL analysis but is not installed")
      res.list[["rhdf5"]] <- "rhdf5"
    }

    if (!is.element("HDF5Array", installed.packages()[, 1])) {
      message(" -- HDF5Array is required for SPATIAL analysis but is not installed")
      res.list[["HDF5Array"]] <- "HDF5Array"
    }
  }

  ### Final message

  if (length(res.list) > 0) {
    message(paste0("               "))
    message(paste0("Check out ", "'https://immunedynamics.github.io/spectre/getting-started/'", " for help with installation"))
  }

  if (length(res.list) == 0) {
    message(" -- All packages successfully installed.")
    message(paste0("               "))
    message(paste0("Check out ", "'https://immunedynamics.github.io/spectre/'", " for protocols"))
  }
}
