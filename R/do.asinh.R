#' Run ArcSinh transformation
#'
#' Transform data in selected columns using ArcSinh transformation with a specified co-factor.
#'
#' @seealso \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' @param dat NO DEFAULT. data.table input sample.
#' @param use.cols NO DEFAULT. Vector of character column names -- these columns will be transformed and added to the data.table as new columns.
#' @param cofactor DEFAULT = 5. Co-factor to use for arcsinh transformation. Can be string of co-factors that align with columns in 'use.cols', just make sure each column is assigned a co-factor (repeat values if necessary).
#' @param auto.cofactor DEFAULT = FALSE. Options are 'flowVS' or 'top10'. Automatically calculates co-factor for each column. 'flowVS' uses the flowVS::estParamFlowVS (takes time). 'top10' uses the 90th percentile of values that are negative (based on absolute value) as the co-factor, and if values are all positive it will use the 10th percentile as co-factor. If 'auto.cofactor' not FALSE, it will override the 'cofactor' argument and force 'append.cf' = TRUE.
#' @param append.cf DEFAULT = FALSE. Appends the co-factor used to the end of the name of the transformed columns.
#' @param reduce.noise DEFAULT = FALSE. This is an experimental calculation which should reduce noise from negative values. Use with caution.
#' @param digits DEFAULT = NULL. Number of decimal places as a limit, not used if NULL. Values beyond will be rounded. Equal to the number or less (i.e. if 9 is used, but only 5 digits are present, then 5 digits will be used). Important to control for small floating point error differences in different OS.
#'
#' @return A data.table with new columns added, that contain the asinh transformed data.
#'
#' @usage do.asinh(dat, use.cols)
#'
#' @import data.table
#'
#' @export do.asinh

do.asinh <- function(dat,
                     use.cols,
                     cofactor = 5,
                     auto.cofactor = FALSE, 
                     append.cf = FALSE,
                     reduce.noise = FALSE,
                     digits = NULL) {
  
  ### Force 'append.cf' = TRUE if using auto.cofactor
  if (auto.cofactor != FALSE) {
    append.cf <- TRUE
  }

  ### Setup data
  value <- dat[, use.cols, with = FALSE]

  ### Numeric checks
  if (isFALSE(all(sapply(value, is.numeric)))) {
    message("It appears that one column in your dataset is non numeric")
    print(sapply(value, is.numeric))
    stop("do.asinh stopped")
  }

  ### Optional noise reduction
  # https://github.com/JinmiaoChenLab/cytofkit/issues/71
  if (reduce.noise == TRUE) {
    message("This noise reduction function is experimental, and should be used with caution")
    value <- value - 1
    loID <- which(value < 0)
    if (length(loID) > 0) {
      value[loID] <- rnorm(length(loID), mean = 0, sd = 0.01)
    }
  }

  ### Arcsinh calculation
  if (auto.cofactor == FALSE) {
    if(length(cofactor) == 1) {
      value <- as.matrix(value)
      value <- value / cofactor
      value <- asinh(value) # value <- log(value + sqrt(value^2 + 1))
      
      if(!is.null(digits)){
        value <- round(value, digits = digits)
      }
      
      value <- data.table::as.data.table(value)
    } else {
      # Check number of co-factors is equal to number of columns
      if (length(cofactor) != length(use.cols)) {
        message("Number of cofactors is not equal to number of columns. Repeat cofactor values if needed or only input one cofactor for all columns.")
        stop("do.asinh stopped")
      }
      
      value <- as.matrix(value)
      
      # Divide each column by the relevant optimised cofactor
      #https://stackoverflow.com/questions/48151278/dividing-columns-of-a-matrix-by-elements-of-a-vector
      value <- sweep(value, 2, cofactor, FUN = '/')
      
      # Calculate asinh
      value <- asinh(value)
      
      value <- data.table::as.data.table(value)
      
      if(!is.null(digits)){
        value <- round(value, digits = digits)
    }
    }
  }
  
  ### Automatic cofactor calculation using 'flowVS'
  if (auto.cofactor == 'flowVS') {
    ## Create flowFrame metadata (column names with descriptions) plus flowFrame
    metadata <- data.frame(name = dimnames(value)[[2]], desc = paste("column", dimnames(value)[[2]], "from dataset"))
    dat.ff <- new("flowFrame",
                  exprs = as.matrix(value), # in order to create a flow frame, data needs to be read as matrix
                  parameters = Biobase::AnnotatedDataFrame(metadata)
    )
    
    # Needs to be a 'flowSet' object...
    dat.ff <- flowCore::flowSet(dat.ff)
    
    ## Calculate cofactors
    cofactors <- flowVS::estParamFlowVS(dat.ff, channels = use.cols)
    cofactors
    
    ## Do arcsinh transformation
    # Prepare data
    value <- dat[, use.cols, with = FALSE]
    value <- as.matrix(value)
    
    # Divide each column by the relevant optimised cofactor
    #https://stackoverflow.com/questions/48151278/dividing-columns-of-a-matrix-by-elements-of-a-vector
    value <- sweep(value, 2, cofactors, FUN = '/')
    
    # Calculate asinh
    value <- asinh(value)
    
    value <- data.table::as.data.table(value)
    
    if(!is.null(digits)){
      value <- round(value, digits = digits)
    }
    
  }
    
  ### Automatic cofactor calculation using 'top10'
    if (auto.cofactor == 'top10') {
      value <- dat[, use.cols, with = FALSE]
      
      ## Calculate cofactors
      cofactors.list <- lapply(value, function(x) {
        if(min(x)<0) {
          abs.dat <- abs(x[x<0])
          
          quantile(abs.dat,
                   probs = 0.9,
                   na.rm = TRUE)
          
        } else {
          # median(x)
          quantile(x,
                   probs = 0.1,
                   na.rm = TRUE)
        }
      })
        
        cofactors.list <- unlist(cofactors.list)
        cofactors.list
        
        ## Do arcsinh transformation
        # Prepare data
        value <- dat[, use.cols, with = FALSE]
        value <- as.matrix(value)
        
        # Divide each column by the relevant optimised cofactor
        #https://stackoverflow.com/questions/48151278/dividing-columns-of-a-matrix-by-elements-of-a-vector
        value <- sweep(value, 2, cofactors.list, FUN = '/')
        
        # Calculate asinh
        value <- asinh(value)
        
        value <- data.table::as.data.table(value)
      
      if(!is.null(digits)){
        value <- round(value, digits = digits)
      }
      
    }
  
  ### Options to append the CF used
  if (append.cf == TRUE) {
    if(auto.cofactor == FALSE) {
      if (length(use.cols) > 1) {
        names(value) <- paste0(names(value), "_asinh_cf", cofactor)
      }
      if (length(use.cols) == 1) {
        names(value) <- paste0(use.cols, "_asinh_cf", cofactor)
      }
    }
    
    if(auto.cofactor == 'flowVS') {
      if (length(use.cols) > 1) {
        names(value) <- paste0(names(value), "_asinh_cf", round(cofactors))
      }
      if (length(use.cols) == 1) {
        names(value) <- paste0(use.cols, "_asinh_cf", round(cofactors))
      }
    }
    
    if(auto.cofactor == "top10") {
      if (length(use.cols) > 1) {
        names(value) <- paste0(names(value), "_asinh_cf", round(cofactors.list))
      }
      if (length(use.cols) == 1) {
        names(value) <- paste0(use.cols, "_asinh_cf", round(cofactors.list))
      }
    }
  }

  if (append.cf == FALSE) {
    if (length(use.cols) > 1) {
      names(value) <- paste0(names(value), "_asinh")
    }
    if (length(use.cols) == 1) {
      names(value) <- paste0(use.cols, "_asinh")
    }
  }
    
  ### Message regarding auto.cofactor
  if (auto.cofactor != FALSE) {
    message("When using 'auto.cofactor', always check your plots! You can always manually set your co-factors.")
  }

  ### Wrap up
  dat <- cbind(dat, value)
  return(dat)
}
