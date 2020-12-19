#' Using logicle to transform data
#' 
#' The function transform data using logicle transformation.
#' It makes use of flowCore's implementation of logicle transformation, and can
#' either automatically infer the transformation function based on your data,
#' or calculate it based on what you want by looking at the default or overriden \code{linearisation.width, max.scale.val, full.transform.width, additional.negative.range}. 
#' Description of some of the parameter is adapted from flowCore's vignette.
#' For more information on what logicle transformation does, please read the manuscript in the references section.
#' 
#' @references
#' Parks, David R., Mario Roederer, and Wayne A. Moore. "A new “Logicle” display method avoids deceptive effects of logarithmic scaling for low signals and compensated data." Cytometry Part A: The Journal of the International Society for Analytical Cytology 69.6 (2006): 541-551. 
#' 
#' @param dat NO DEFAULT. data.table. Input data
#' @param use.cols NO DEFAULT. Vector of column names to transform
#' @param auto.infer.function Default = TRUE. 
#'   Automatically infer the logical transformation function based on your data.
#'   If this is set to FALSE, default or overriden 
#'   \code{linearisation.width, max.scale.val, full.transform.width, additional.negative.range}
#'   will be used instead to calculate the logical transformation function.
#' @param linearisation.width Default = 1.2. Linearisation width in asymptotic decades.
#'   This must be > 0 and determines the slope of transformation at zero.
#'   It can be estimated using equation:
#'   \eqn{(m-log10(max.scale.val/|r|))/2}
#'   where r is the most negative value to included in the display.
#'   See \code{max.scale.val} below.
#' @param max.scale.val Default = 262144. Maximum scale data value. 
#'   It can be 10,000 for common 4 decade data or 262144 (the default value) for 18 bit data range.
#'   This must be greater than 0
#' @param full.transform.width Default = 4.5. The full width of the transformed values in asymptotic decades
#'   This must be greater than 0
#' @param additional.negative.range Default = 0. Additional negative range to be included in 
#'   the transformed value in asymptotic decades.
#'   Value greater than 0 will bring the prescribed additional range into the transformed values
#'   
#' @usage 
#' do.logicle(dat, use.cols, linearisation.width, max.scale.val, full.transform.width, additional.negative.range)
#'   
#' @examples 
#' library(Spectre)
#' dat = demo.asinh
#' use.cols = names(dat)[c(2:10)]
#' dat_lgcl = do.logicle(dat, use.cols = use.cols)
#' 
#' @author
#' Givanna Putri, \email{givanna.haryonoputri@@sydney.edu.au}
#'   
#' @export

do.logicle <- function(dat, 
                       use.cols,
                       auto.infer.function = TRUE,
                       linearisation.width = 1.2,
                       max.scale.val = 262144,
                       full.transform.width = 4.5,
                       additional.negative.range = 0) {
  
  require(flowCore)
  require(Biobase)
  require(data.table)
  
  # check selected columns are numeric
  if(any(unlist(lapply(dat[, use.cols, with = FALSE], is.numeric)) == FALSE)) {
    stop('Non-numeric columns are selected for transformation. Please exclude them.')
  }
  
  values <- dat[,use.cols, with=FALSE]
  
  val.mat <- as.matrix(values)
  
  ff <- new("flowFrame", exprs = val.mat)
  
  if (auto.infer.function) {
    # Automatic estimated logicle function
    message("Automatically estimating the logicle transformation function based on input data")
    lgcl <- estimateLogicle(ff, channels = use.cols) 
  } else {
    
    # check parameter conditions, make sure they are greater than 0
    if (is.null(max.scale.val) || max.scale.val <= 0) {
      stop('max.scale.val must be greater than 0')
    }
    if (is.null(linearisation.width) || linearisation.width <= 0) {
      stop('linearisation.width must be greater than 0')
    }
    if (is.null(full.transform.width) || full.transform.width <= 0) {
      stop('full.transform.width must be greater than 0')
    }
    
    # User defined logicle function
    message("Formulating logicle transformation function")
    trans_lgcl <- logicleTransform(w = linearisation.width, 
                                   t = max.scale.val, 
                                   m = full.transform.width, 
                                   a = additional.negative.range)
    lgcl <- transformList(use.cols, trans_lgcl)
  }
  
  # Do transformation
  message("Transforming data")
  trans.val <- flowCore::transform(ff, lgcl)
  
  message("Converting data back to data.table")
  trans.val.dt <- data.table(exprs(trans.val))
  names(trans.val.dt) <- paste(names(trans.val.dt), "logicle", sep = "_")
  
  dat_copy <- cbind(dat, trans.val.dt)
  return(dat_copy)
}