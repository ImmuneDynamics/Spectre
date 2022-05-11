#' package.install - a function to install packages required for Spectre.
#'
#' This function allows you to check to see if all the common use packages dependencies for Spectre are installed. Will only install if the package has not been installed, will not update packages.
#'
#' @return returns an error message if one of the common use packages are not installed. Proceeds in order of package importance, and only the first error message encountered will be returned.
#'
#' @param type DEFAULT = "general". If "general", then checks for the packages required for general Spectre usage. If "spatial", then checks for additional packages required for spatial analysis.
#' @param update DEFAULT = FALSE. If FALSE, will only install packages that are not already installed -- no updates of packages will be performed. If TRUE, will install and update packages.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' @seealso See \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
#'
#' @usage package.install()
#'
#' @examples
#' package.install()
#' package.install(type = "spatial")
#'
#' @export

package.install <- function(type = "general", update = FALSE) {

  ## Update = FALSE
  if (isFALSE(update)) {
    if (!require("data.table")) {
      install.packages("data.table")
    }
    if (!require("plyr")) {
      install.packages("plyr")
    }
    if (!require("dplyr")) {
      install.packages("dplyr")
    }
    if (!require("tidyr")) {
      install.packages("tidyr")
    }
    if (!require("rstudioapi")) {
      install.packages("rstudioapi")
    }
    if (!require("Rtsne")) {
      install.packages("Rtsne")
    }
    if (!require("umap")) {
      install.packages("umap")
    }
    if (!require("ggplot2")) {
      install.packages("ggplot2")
    }
    if (!require("ggthemes")) {
      install.packages("ggthemes")
    }
    if (!require("scales")) {
      install.packages("scales")
    }
    if (!require("colorRamps")) {
      install.packages("colorRamps")
    }
    if (!require("RColorBrewer")) {
      install.packages("")
    }
    if (!require("gridExtra")) {
      install.packages("gridExtra")
    }
    if (!require("ggpointdensity")) {
      install.packages("ggpointdensity")
    }
    if (!require("pheatmap")) {
      install.packages("pheatmap")
    }
    if (!require("ggpubr")) {
      install.packages("ggpubr")
    }
    if (!require("factoextra")) {
      install.packages("factoextra")
    }
    if (!require("reticulate")) {
      install.packages("reticulate")
    }
    if (!require("caret")) {
      install.packages("caret")
    }
    if (!require("class")) {
      install.packages("class")
    }

    ## Install BiocManager to download packages from Bioconductor
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }

    ## Download additional BioConductor packages
    if (!require("flowCore")) {
      BiocManager::install("flowCore")
    }
    if (!require("Biobase")) {
      BiocManager::install("Biobase")
    }
    if (!require("flowViz")) {
      BiocManager::install("flowViz")
    }
    if (!require("FlowSOM")) {
      BiocManager::install("FlowSOM")
    }

    if (type == "spatial") {
      if (!require("raster")) {
        install.packages("raster")
      }
      if (!require("tiff")) {
        install.packages("tiff")
      }
      if (!require("rgeos")) {
        install.packages("rgeos")
      }
      if (!require("exactextractr")) {
        install.packages("exactextractr")
      }
      if (!require("sp")) {
        install.packages("sp")
      }
      if (!require("sf")) {
        install.packages("sf")
      }
      if (!require("stars")) {
        install.packages("stars")
      }
      if (!require("qs")) {
        install.packages("qs")
      }
      if (!require("s2")) {
        install.packages("s2")
      }
      if (!require("qs")) {
        install.packages("HDF5Array")
      }
      if (!require("s2")) {
        install.packages("hdf5r")
      }
    }
  }

  ## Update = TRUE
  if (isTRUE(update)) {
    install.packages("data.table")
    install.packages("plyr")
    install.packages("dplyr")
    install.packages("tidyr")
    install.packages("rstudioapi")
    install.packages("Rtsne")
    install.packages("umap")
    install.packages("ggplot2")
    install.packages("ggthemes")
    install.packages("scales")
    install.packages("colorRamps")
    install.packages("RColorBrewer")
    install.packages("gridExtra")
    install.packages("ggpointdensity")
    install.packages("pheatmap")
    install.packages("ggpubr")
    install.packages("factoextra")
    install.packages("reticulate")
    install.packages("caret")
    install.packages("class")

    ## Install BiocManager to download packages from Bioconductor
    install.packages("BiocManager")

    ## Download additional BioConductor packages
    BiocManager::install("flowCore")
    BiocManager::install("Biobase")
    BiocManager::install("flowViz")
    BiocManager::install("FlowSOM")

    if (type == "spatial") {
      install.packages("raster")
      install.packages("tiff")
      install.packages("rgeos")
      install.packages("exactextractr")
      install.packages("sp")
      install.packages("sf")
      install.packages("stars")
      install.packages("qs")
      install.packages("s2")
      install.packages("HDF5Array")
      install.packages("rhdf5")
    }
  }
}
