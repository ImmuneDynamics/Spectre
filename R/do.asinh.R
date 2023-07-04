#' Run ArcSinh transformation
#'
#' Transform data in selected columns using ArcSinh transformation with a specified co-factor.
#'
#' @seealso \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' @param dat NO DEFAULT. data.table input sample.
#' @param use.cols NO DEFAULT. Vector of character column names -- these columns will be transformed and added to the data.table as new columns.
#' @param cofactor DEFAULT = 5. Co-factor to use for arcsinh transformation.
#' @param append.cf DEFAULT = FALSE. Appends the co-factor used to the end of the name of the transformed columns.
#' @param reduce.noise DEFAULT = FALSE. This is an experimental calculation which should reduce noise from negative values. Use with caution.
#' @param digits DEFAULT = 9. Number of decimal places as a limit. Values beyond will be rounded. Equal to the number or less (i.e. if 9 is used, but only 5 digits are present, then 5 digits will be used). Important to control for small floating point error differences in different OS.
#'
#' @return A data.table with new columns added, that contain the asinh transformed data.
#'
#' @usage do.asinh(dat, use.cols)
#'
#' @import data.table
#'
#' @export do.asinh

do.asinh <- function(dat,
                     use.cols = NULL,
                     cofactor = 5,
                     assay = NULL,
                     clip.min = NULL,
                     clip.max = NULL,
                     name = '_asinh',
                     digits = 9) {
  
  ### Evaluation
  
      if(class(dat)[1] == 'data.table'){
        if(!is.null(use.cols)){
          value <- as.matrix(dat)
        } else {
          value <- as.matrix(dat[,use.cols, with = FALSE])
        }
      }
  
      if(class(dat)[1] == 'spectre'){
        if(is.null(assay)){
          stop('assay required')
        }
        value <- as.matrix(dat@data[[assay]])
        if(!is.null(use.cols)){
          value <- value[,use.cols]
        }
      }
  
      if(class(dat)[1] == 'Seurat'){
        if(is.null(assay)){
          stop('assay required')
        }
        value <- as.matrix(dat@assays[[assay]]@counts)
        if(!is.null(use.cols)){
          value <- value[,use.cols]
          value <- t(value)
        }
      }

  ### Numeric checks
  
      if (isFALSE(all(sapply(value, is.numeric)))) {
        message("It appears that one column in your dataset is non numeric")
        print(sapply(value, is.numeric))
        stop("do.asinh stopped")
      }
  
  ### Arcsinh calculation
  
      value <- value / cofactor
      value <- asinh(value) # value <- log(value + sqrt(value^2 + 1))
      value <- round(value, digits = digits)

  ### Clipping
      
      if(!is.null(clip.min)){
        
      }
      
      if(!is.null(clip.max)){
        
      }
      
  ### Wrap up
      
      if(class(dat)[1] == 'data.table'){
        names(value) <- paste0(names(value), '_', name)
        # Add test to see if those cols are there, if so, replace
        dat <- rbind(dat, value)
      }
      if(class(dat)[1] == 'spectre'){
        dat@data[[paste0(assay, '_', name)]] <- value
      }
      if(class(dat)[1] == 'Seurat'){
        dat@assays[[assay]]@data <- value
      }
      
  ### Return

      return(dat)

}
