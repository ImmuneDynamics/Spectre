#' Generic funcion to do ArcSinh transformation
#'
#' Transform data in selected columns using ArcSinh transformation with a specified co-factor.
#'
#' @seealso \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' @param dat NO DEFAULT. 
#' Either a data.table or a SpectreObject.
#' @param use.cols NO DEFAULT. 
#' A vector of character column names.
#' These columns will be transformed and added to the data.table as new columns.
#' If `dat` is a SpectreObject, the values in this parameter must exists in the
#' columns of the slot specified in `slot` parameter!
#' @param cofactor DEFAULT = 5. Co-factor to use for arcsinh transformation.
#' @param append.cf DEFAULT = FALSE. Appends the co-factor used to the end of 
#' the name of the transformed columns.
#' @param reduce.noise DEFAULT = FALSE. This is an experimental calculation 
#' which should reduce noise from negative values. Use with caution.
#' @param digits DEFAULT = NULL. Number of decimal places as a limit, not used if NULL. 
#' Values beyond will be rounded. 
#' Equal to the number or less (i.e. if 9 is used, but only 5 digits are present, 
#' then 5 digits will be used). Important to control for small floating point 
#' error differences in different OS.
#' @param slot_name DEFAULT = NULL.
#' A mandatory parameter for SpectreObject method.
#' Otherwise not required.
#' It specifies which slot to apply the arcsinh transformation to.
#' Can be a vector, which then will apply the arcsinh to multiple slots
#' or a character, which will then apply the arcsinh only to that slot.
#' If left as it is (NULL), no transformation will be performed!
#' See details why this is set by default to NULL.
#' @param verbose DEFAULT = TRUE
#' If TRUE, the function will print progress updates as it executes.
#' 
#' @details
#' slot_name parameter is set to NULL by default because if the user accidentally
#' forgot to pass the slot name to apply the transformation to, R gives gibberish error. 
#' It'll make it so hard to debug.
#' By setting it to NULL by default, if the user forgot to pass the slot name,
#' we can throw out a warning message saying nothing was done because no slot
#' was specified.
#' It also overall makes it easier to debug the code.
#' Could use try and catch. In the future perhaps.
#' 
#'
#' @return A data.table with new columns added, that contain the asinh transformed data.
#'
#' @usage do.asinh(dat, use.cols)
#'
#' @export
setGeneric("do.asinh", function(dat,
                                use.cols,
                                slot_name = NULL,
                                cofactor = 5,
                                append.cf = FALSE,
                                reduce.noise = FALSE,
                                digits = NULL,
                                verbose = TRUE) {
    standardGeneric("do.asinh")
    
})