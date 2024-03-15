#' Do ArcSinh transformation
#'
#' Transform data in selected columns using ArcSinh transformation with a specified co-factor.
#'
#' @seealso \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' @param dat NO DEFAULT. 
#' Either a data.table or a Spectre object to apply arc-sinh transformation to.
#' @param use.cols NO DEFAULT. 
#' A vector of character column names to apply arc-sinh transformation to.
#' @param cofactor DEFAULT = 5. Co-factor to use for arcsinh transformation.
#' Supply a data.table with column marker and cofactor if you want to apply different co-factor to different markers.
#' See examples.
#' @param append.cf DEFAULT = FALSE. Appends the co-factor used to the end of 
#' the name of the transformed columns.
#' @param reduce.noise DEFAULT = FALSE. This is an experimental calculation 
#' which should reduce noise from negative values. Use with caution.
#' @param digits DEFAULT = NULL. Number of decimal places as a limit, not used if NULL. 
#' Values beyond will be rounded. 
#' Equal to the number or less (i.e. if 9 is used, but only 5 digits are present, 
#' then 5 digits will be used). Important to control for small floating point 
#' error differences in different OS.
#' @param verbose DEFAULT = TRUE.
#' If TRUE, the function will print progress updates as it executes.
#' 
#' @examples
#' library(data.table)
#' # Assuming dat is a data.table
#' dat <- data.table(NK11=rnorm(10, 2), CD3=rnorm(10,1))
#' # Default co-factor 5 will be used for both NK11 and CD3
#' dat_asinh <- do.asinh(dat, use.cols=c("NK11", "CD3"))
#' dat_asinh
#' 
#' # Apply asinh to only NK11
#' dat_asinh <- do.asinh(dat, use.cols=c("NK11"))
#' dat_asinh
#' 
#' # Apply different co-factor to the markers
#' cofactor_dt <- data.table(marker=c("NK11", "CD3"), cofactor=c(5,10))
#' cofactor_dt
#' dat_asinh <- do.asinh(dat, use.cols=c("NK11", "CD3"), cofactor=cofactor_dt)
#' dat_asinh
#' 
#' 
#' @export
setGeneric("do.asinh", function(dat,
                                use.cols,
                                cofactor = 5,
                                append.cf = FALSE,
                                reduce.noise = FALSE,
                                digits = NULL,
                                verbose = TRUE,
                                ...) {
    standardGeneric("do.asinh")
    
})