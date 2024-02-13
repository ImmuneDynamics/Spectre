#' do.filter 
#' 
#' This function allows filtering of a data.table using multiple match values -- all cells that contain any of the match values will be filtered, and provided in a new data.table.
#'
#' @param dat NO DEFAULT. A data.table
#' @param use.col DEFAULT = NULL. The column to use for re-ordering
#' @param values DEFAULT = NULL. A vector of values to use for filtering -- all cells that contain any of the match values will be filtered, and provided in a new data.table. Up to 35 match values can be provided
#'
#' @usage do.filter(dat, use.col, values)
#'
#' @author
#' Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://github.com/ImmuneDynamics/Spectre}.
#'
#'
#' @export

  do.filter <- function(dat,
                        use.col,
                        values) {
    # require: data.table
    
    # TODO: can be improved i think.
  
    ### Demo data
  
    # dat <- Spectre::demo.clustered
    # use.col <- 'Population'
    #
    # unique(dat[[use.col]])
    #
    # values <- c("Microglia", "Infil Macrophages")
  
    ### Setup
  
    n.values <- length(values)
  
    ### Warnings
  
    if (n.values > 35) {
      stop("do.filter function cannot take more than 35 target values")
    }
  
    ### Filtering
  
    if (n.values == 1) {
      res <- dat[dat[[use.col]] == values[1], ]
    }
  
    if (n.values == 2) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2], ]
    }
  
    if (n.values == 3) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2] |
        dat[[use.col]] == values[3], ]
    }
  
    if (n.values == 4) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2] |
        dat[[use.col]] == values[3] |
        dat[[use.col]] == values[4], ]
    }
  
    if (n.values == 5) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2] |
        dat[[use.col]] == values[3] |
        dat[[use.col]] == values[4] |
        dat[[use.col]] == values[5], ]
    }
  
    if (n.values == 6) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2] |
        dat[[use.col]] == values[3] |
        dat[[use.col]] == values[4] |
        dat[[use.col]] == values[5] |
        dat[[use.col]] == values[6], ]
    }
  
    if (n.values == 7) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2] |
        dat[[use.col]] == values[3] |
        dat[[use.col]] == values[4] |
        dat[[use.col]] == values[5] |
        dat[[use.col]] == values[6] |
        dat[[use.col]] == values[7], ]
    }
  
    if (n.values == 8) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2] |
        dat[[use.col]] == values[3] |
        dat[[use.col]] == values[4] |
        dat[[use.col]] == values[5] |
        dat[[use.col]] == values[6] |
        dat[[use.col]] == values[7] |
        dat[[use.col]] == values[8], ]
    }
  
    if (n.values == 9) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2] |
        dat[[use.col]] == values[3] |
        dat[[use.col]] == values[4] |
        dat[[use.col]] == values[5] |
        dat[[use.col]] == values[6] |
        dat[[use.col]] == values[7] |
        dat[[use.col]] == values[8] |
        dat[[use.col]] == values[9], ]
    }
  
    if (n.values == 10) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2] |
        dat[[use.col]] == values[3] |
        dat[[use.col]] == values[4] |
        dat[[use.col]] == values[5] |
        dat[[use.col]] == values[6] |
        dat[[use.col]] == values[7] |
        dat[[use.col]] == values[8] |
        dat[[use.col]] == values[9] |
        dat[[use.col]] == values[10], ]
    }
  
    if (n.values == 11) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2] |
        dat[[use.col]] == values[3] |
        dat[[use.col]] == values[4] |
        dat[[use.col]] == values[5] |
        dat[[use.col]] == values[6] |
        dat[[use.col]] == values[7] |
        dat[[use.col]] == values[8] |
        dat[[use.col]] == values[9] |
        dat[[use.col]] == values[10] |
        dat[[use.col]] == values[11], ]
    }
  
    if (n.values == 12) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2] |
        dat[[use.col]] == values[3] |
        dat[[use.col]] == values[4] |
        dat[[use.col]] == values[5] |
        dat[[use.col]] == values[6] |
        dat[[use.col]] == values[7] |
        dat[[use.col]] == values[8] |
        dat[[use.col]] == values[9] |
        dat[[use.col]] == values[10] |
        dat[[use.col]] == values[11] |
        dat[[use.col]] == values[12], ]
    }
  
    if (n.values == 13) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2] |
        dat[[use.col]] == values[3] |
        dat[[use.col]] == values[4] |
        dat[[use.col]] == values[5] |
        dat[[use.col]] == values[6] |
        dat[[use.col]] == values[7] |
        dat[[use.col]] == values[8] |
        dat[[use.col]] == values[9] |
        dat[[use.col]] == values[10] |
        dat[[use.col]] == values[11] |
        dat[[use.col]] == values[12] |
        dat[[use.col]] == values[13], ]
    }
  
    if (n.values == 14) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2] |
        dat[[use.col]] == values[3] |
        dat[[use.col]] == values[4] |
        dat[[use.col]] == values[5] |
        dat[[use.col]] == values[6] |
        dat[[use.col]] == values[7] |
        dat[[use.col]] == values[8] |
        dat[[use.col]] == values[9] |
        dat[[use.col]] == values[10] |
        dat[[use.col]] == values[11] |
        dat[[use.col]] == values[12] |
        dat[[use.col]] == values[13] |
        dat[[use.col]] == values[14], ]
    }
  
    if (n.values == 15) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2] |
        dat[[use.col]] == values[3] |
        dat[[use.col]] == values[4] |
        dat[[use.col]] == values[5] |
        dat[[use.col]] == values[6] |
        dat[[use.col]] == values[7] |
        dat[[use.col]] == values[8] |
        dat[[use.col]] == values[9] |
        dat[[use.col]] == values[10] |
        dat[[use.col]] == values[11] |
        dat[[use.col]] == values[12] |
        dat[[use.col]] == values[13] |
        dat[[use.col]] == values[14] |
        dat[[use.col]] == values[15], ]
    }
  
    if (n.values == 16) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2] |
        dat[[use.col]] == values[3] |
        dat[[use.col]] == values[4] |
        dat[[use.col]] == values[5] |
        dat[[use.col]] == values[6] |
        dat[[use.col]] == values[7] |
        dat[[use.col]] == values[8] |
        dat[[use.col]] == values[9] |
        dat[[use.col]] == values[10] |
        dat[[use.col]] == values[11] |
        dat[[use.col]] == values[12] |
        dat[[use.col]] == values[13] |
        dat[[use.col]] == values[14] |
        dat[[use.col]] == values[15] |
        dat[[use.col]] == values[16], ]
    }
  
    if (n.values == 17) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2] |
        dat[[use.col]] == values[3] |
        dat[[use.col]] == values[4] |
        dat[[use.col]] == values[5] |
        dat[[use.col]] == values[6] |
        dat[[use.col]] == values[7] |
        dat[[use.col]] == values[8] |
        dat[[use.col]] == values[9] |
        dat[[use.col]] == values[10] |
        dat[[use.col]] == values[11] |
        dat[[use.col]] == values[12] |
        dat[[use.col]] == values[13] |
        dat[[use.col]] == values[14] |
        dat[[use.col]] == values[15] |
        dat[[use.col]] == values[16] |
        dat[[use.col]] == values[17], ]
    }
  
    if (n.values == 18) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2] |
        dat[[use.col]] == values[3] |
        dat[[use.col]] == values[4] |
        dat[[use.col]] == values[5] |
        dat[[use.col]] == values[6] |
        dat[[use.col]] == values[7] |
        dat[[use.col]] == values[8] |
        dat[[use.col]] == values[9] |
        dat[[use.col]] == values[10] |
        dat[[use.col]] == values[11] |
        dat[[use.col]] == values[12] |
        dat[[use.col]] == values[13] |
        dat[[use.col]] == values[14] |
        dat[[use.col]] == values[15] |
        dat[[use.col]] == values[16] |
        dat[[use.col]] == values[17] |
        dat[[use.col]] == values[18], ]
    }
  
    if (n.values == 19) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2] |
        dat[[use.col]] == values[3] |
        dat[[use.col]] == values[4] |
        dat[[use.col]] == values[5] |
        dat[[use.col]] == values[6] |
        dat[[use.col]] == values[7] |
        dat[[use.col]] == values[8] |
        dat[[use.col]] == values[9] |
        dat[[use.col]] == values[10] |
        dat[[use.col]] == values[11] |
        dat[[use.col]] == values[12] |
        dat[[use.col]] == values[13] |
        dat[[use.col]] == values[14] |
        dat[[use.col]] == values[15] |
        dat[[use.col]] == values[16] |
        dat[[use.col]] == values[17] |
        dat[[use.col]] == values[18] |
        dat[[use.col]] == values[19], ]
    }
  
    if (n.values == 20) {
      res <- dat[dat[[use.col]] == values[1] |
        dat[[use.col]] == values[2] |
        dat[[use.col]] == values[3] |
        dat[[use.col]] == values[4] |
        dat[[use.col]] == values[5] |
        dat[[use.col]] == values[6] |
        dat[[use.col]] == values[7] |
        dat[[use.col]] == values[8] |
        dat[[use.col]] == values[9] |
        dat[[use.col]] == values[10] |
        dat[[use.col]] == values[11] |
        dat[[use.col]] == values[12] |
        dat[[use.col]] == values[13] |
        dat[[use.col]] == values[14] |
        dat[[use.col]] == values[15] |
        dat[[use.col]] == values[16] |
        dat[[use.col]] == values[17] |
        dat[[use.col]] == values[18] |
        dat[[use.col]] == values[19] |
        dat[[use.col]] == values[20], ]
    }
    
    if (n.values == 21) {
      res <- dat[dat[[use.col]] == values[1] |
                   dat[[use.col]] == values[2] |
                   dat[[use.col]] == values[3] |
                   dat[[use.col]] == values[4] |
                   dat[[use.col]] == values[5] |
                   dat[[use.col]] == values[6] |
                   dat[[use.col]] == values[7] |
                   dat[[use.col]] == values[8] |
                   dat[[use.col]] == values[9] |
                   dat[[use.col]] == values[10] |
                   dat[[use.col]] == values[11] |
                   dat[[use.col]] == values[12] |
                   dat[[use.col]] == values[13] |
                   dat[[use.col]] == values[14] |
                   dat[[use.col]] == values[15] |
                   dat[[use.col]] == values[16] |
                   dat[[use.col]] == values[17] |
                   dat[[use.col]] == values[18] |
                   dat[[use.col]] == values[19] |
                   dat[[use.col]] == values[20] |
                   dat[[use.col]] == values[21], ]
    }
    
    if (n.values == 22) {
      res <- dat[dat[[use.col]] == values[1] |
                   dat[[use.col]] == values[2] |
                   dat[[use.col]] == values[3] |
                   dat[[use.col]] == values[4] |
                   dat[[use.col]] == values[5] |
                   dat[[use.col]] == values[6] |
                   dat[[use.col]] == values[7] |
                   dat[[use.col]] == values[8] |
                   dat[[use.col]] == values[9] |
                   dat[[use.col]] == values[10] |
                   dat[[use.col]] == values[11] |
                   dat[[use.col]] == values[12] |
                   dat[[use.col]] == values[13] |
                   dat[[use.col]] == values[14] |
                   dat[[use.col]] == values[15] |
                   dat[[use.col]] == values[16] |
                   dat[[use.col]] == values[17] |
                   dat[[use.col]] == values[18] |
                   dat[[use.col]] == values[19] |
                   dat[[use.col]] == values[20] |
                   dat[[use.col]] == values[21] |
                   dat[[use.col]] == values[22], ]
    }
    
    if (n.values == 23) {
      res <- dat[dat[[use.col]] == values[1] |
                   dat[[use.col]] == values[2] |
                   dat[[use.col]] == values[3] |
                   dat[[use.col]] == values[4] |
                   dat[[use.col]] == values[5] |
                   dat[[use.col]] == values[6] |
                   dat[[use.col]] == values[7] |
                   dat[[use.col]] == values[8] |
                   dat[[use.col]] == values[9] |
                   dat[[use.col]] == values[10] |
                   dat[[use.col]] == values[11] |
                   dat[[use.col]] == values[12] |
                   dat[[use.col]] == values[13] |
                   dat[[use.col]] == values[14] |
                   dat[[use.col]] == values[15] |
                   dat[[use.col]] == values[16] |
                   dat[[use.col]] == values[17] |
                   dat[[use.col]] == values[18] |
                   dat[[use.col]] == values[19] |
                   dat[[use.col]] == values[20] |
                   dat[[use.col]] == values[21] |
                   dat[[use.col]] == values[22] |
                   dat[[use.col]] == values[23], ]
    }
    
    if (n.values == 24) {
      res <- dat[dat[[use.col]] == values[1] |
                   dat[[use.col]] == values[2] |
                   dat[[use.col]] == values[3] |
                   dat[[use.col]] == values[4] |
                   dat[[use.col]] == values[5] |
                   dat[[use.col]] == values[6] |
                   dat[[use.col]] == values[7] |
                   dat[[use.col]] == values[8] |
                   dat[[use.col]] == values[9] |
                   dat[[use.col]] == values[10] |
                   dat[[use.col]] == values[11] |
                   dat[[use.col]] == values[12] |
                   dat[[use.col]] == values[13] |
                   dat[[use.col]] == values[14] |
                   dat[[use.col]] == values[15] |
                   dat[[use.col]] == values[16] |
                   dat[[use.col]] == values[17] |
                   dat[[use.col]] == values[18] |
                   dat[[use.col]] == values[19] |
                   dat[[use.col]] == values[20] |
                   dat[[use.col]] == values[21] |
                   dat[[use.col]] == values[22] |
                   dat[[use.col]] == values[23] |
                   dat[[use.col]] == values[24], ]
    }
    
    
    if (n.values == 25) {
      res <- dat[dat[[use.col]] == values[1] |
                   dat[[use.col]] == values[2] |
                   dat[[use.col]] == values[3] |
                   dat[[use.col]] == values[4] |
                   dat[[use.col]] == values[5] |
                   dat[[use.col]] == values[6] |
                   dat[[use.col]] == values[7] |
                   dat[[use.col]] == values[8] |
                   dat[[use.col]] == values[9] |
                   dat[[use.col]] == values[10] |
                   dat[[use.col]] == values[11] |
                   dat[[use.col]] == values[12] |
                   dat[[use.col]] == values[13] |
                   dat[[use.col]] == values[14] |
                   dat[[use.col]] == values[15] |
                   dat[[use.col]] == values[16] |
                   dat[[use.col]] == values[17] |
                   dat[[use.col]] == values[18] |
                   dat[[use.col]] == values[19] |
                   dat[[use.col]] == values[20] |
                   dat[[use.col]] == values[21] |
                   dat[[use.col]] == values[22] |
                   dat[[use.col]] == values[23] |
                   dat[[use.col]] == values[24] |
                   dat[[use.col]] == values[25], ]
    }
    
    if (n.values == 26) {
      res <- dat[dat[[use.col]] == values[1] |
                   dat[[use.col]] == values[2] |
                   dat[[use.col]] == values[3] |
                   dat[[use.col]] == values[4] |
                   dat[[use.col]] == values[5] |
                   dat[[use.col]] == values[6] |
                   dat[[use.col]] == values[7] |
                   dat[[use.col]] == values[8] |
                   dat[[use.col]] == values[9] |
                   dat[[use.col]] == values[10] |
                   dat[[use.col]] == values[11] |
                   dat[[use.col]] == values[12] |
                   dat[[use.col]] == values[13] |
                   dat[[use.col]] == values[14] |
                   dat[[use.col]] == values[15] |
                   dat[[use.col]] == values[16] |
                   dat[[use.col]] == values[17] |
                   dat[[use.col]] == values[18] |
                   dat[[use.col]] == values[19] |
                   dat[[use.col]] == values[20] |
                   dat[[use.col]] == values[21] |
                   dat[[use.col]] == values[22] |
                   dat[[use.col]] == values[23] |
                   dat[[use.col]] == values[24] |
                   dat[[use.col]] == values[25] |
                   dat[[use.col]] == values[26], ]
    }
    
    if (n.values == 27) {
      res <- dat[dat[[use.col]] == values[1] |
                   dat[[use.col]] == values[2] |
                   dat[[use.col]] == values[3] |
                   dat[[use.col]] == values[4] |
                   dat[[use.col]] == values[5] |
                   dat[[use.col]] == values[6] |
                   dat[[use.col]] == values[7] |
                   dat[[use.col]] == values[8] |
                   dat[[use.col]] == values[9] |
                   dat[[use.col]] == values[10] |
                   dat[[use.col]] == values[11] |
                   dat[[use.col]] == values[12] |
                   dat[[use.col]] == values[13] |
                   dat[[use.col]] == values[14] |
                   dat[[use.col]] == values[15] |
                   dat[[use.col]] == values[16] |
                   dat[[use.col]] == values[17] |
                   dat[[use.col]] == values[18] |
                   dat[[use.col]] == values[19] |
                   dat[[use.col]] == values[20] |
                   dat[[use.col]] == values[21] |
                   dat[[use.col]] == values[22] |
                   dat[[use.col]] == values[23] |
                   dat[[use.col]] == values[24] |
                   dat[[use.col]] == values[25] |
                   dat[[use.col]] == values[26] |
                   dat[[use.col]] == values[27], ]
    }
    
    if (n.values == 28) {
      res <- dat[dat[[use.col]] == values[1] |
                   dat[[use.col]] == values[2] |
                   dat[[use.col]] == values[3] |
                   dat[[use.col]] == values[4] |
                   dat[[use.col]] == values[5] |
                   dat[[use.col]] == values[6] |
                   dat[[use.col]] == values[7] |
                   dat[[use.col]] == values[8] |
                   dat[[use.col]] == values[9] |
                   dat[[use.col]] == values[10] |
                   dat[[use.col]] == values[11] |
                   dat[[use.col]] == values[12] |
                   dat[[use.col]] == values[13] |
                   dat[[use.col]] == values[14] |
                   dat[[use.col]] == values[15] |
                   dat[[use.col]] == values[16] |
                   dat[[use.col]] == values[17] |
                   dat[[use.col]] == values[18] |
                   dat[[use.col]] == values[19] |
                   dat[[use.col]] == values[20] |
                   dat[[use.col]] == values[21] |
                   dat[[use.col]] == values[22] |
                   dat[[use.col]] == values[23] |
                   dat[[use.col]] == values[24] |
                   dat[[use.col]] == values[25] |
                   dat[[use.col]] == values[26] |
                   dat[[use.col]] == values[27] |
                   dat[[use.col]] == values[28], ]
    }
    
    if (n.values == 29) {
      res <- dat[dat[[use.col]] == values[1] |
                   dat[[use.col]] == values[2] |
                   dat[[use.col]] == values[3] |
                   dat[[use.col]] == values[4] |
                   dat[[use.col]] == values[5] |
                   dat[[use.col]] == values[6] |
                   dat[[use.col]] == values[7] |
                   dat[[use.col]] == values[8] |
                   dat[[use.col]] == values[9] |
                   dat[[use.col]] == values[10] |
                   dat[[use.col]] == values[11] |
                   dat[[use.col]] == values[12] |
                   dat[[use.col]] == values[13] |
                   dat[[use.col]] == values[14] |
                   dat[[use.col]] == values[15] |
                   dat[[use.col]] == values[16] |
                   dat[[use.col]] == values[17] |
                   dat[[use.col]] == values[18] |
                   dat[[use.col]] == values[19] |
                   dat[[use.col]] == values[20] |
                   dat[[use.col]] == values[21] |
                   dat[[use.col]] == values[22] |
                   dat[[use.col]] == values[23] |
                   dat[[use.col]] == values[24] |
                   dat[[use.col]] == values[25] |
                   dat[[use.col]] == values[26] |
                   dat[[use.col]] == values[27] |
                   dat[[use.col]] == values[28] |
                   dat[[use.col]] == values[29], ]
    }
    
    if (n.values == 30) {
      res <- dat[dat[[use.col]] == values[1] |
                   dat[[use.col]] == values[2] |
                   dat[[use.col]] == values[3] |
                   dat[[use.col]] == values[4] |
                   dat[[use.col]] == values[5] |
                   dat[[use.col]] == values[6] |
                   dat[[use.col]] == values[7] |
                   dat[[use.col]] == values[8] |
                   dat[[use.col]] == values[9] |
                   dat[[use.col]] == values[10] |
                   dat[[use.col]] == values[11] |
                   dat[[use.col]] == values[12] |
                   dat[[use.col]] == values[13] |
                   dat[[use.col]] == values[14] |
                   dat[[use.col]] == values[15] |
                   dat[[use.col]] == values[16] |
                   dat[[use.col]] == values[17] |
                   dat[[use.col]] == values[18] |
                   dat[[use.col]] == values[19] |
                   dat[[use.col]] == values[20] |
                   dat[[use.col]] == values[21] |
                   dat[[use.col]] == values[22] |
                   dat[[use.col]] == values[23] |
                   dat[[use.col]] == values[24] |
                   dat[[use.col]] == values[25] |
                   dat[[use.col]] == values[26] |
                   dat[[use.col]] == values[27] |
                   dat[[use.col]] == values[28] |
                   dat[[use.col]] == values[29] |
                   dat[[use.col]] == values[30], ]
    }
    
    if (n.values == 31) {
      res <- dat[dat[[use.col]] == values[1] |
                   dat[[use.col]] == values[2] |
                   dat[[use.col]] == values[3] |
                   dat[[use.col]] == values[4] |
                   dat[[use.col]] == values[5] |
                   dat[[use.col]] == values[6] |
                   dat[[use.col]] == values[7] |
                   dat[[use.col]] == values[8] |
                   dat[[use.col]] == values[9] |
                   dat[[use.col]] == values[10] |
                   dat[[use.col]] == values[11] |
                   dat[[use.col]] == values[12] |
                   dat[[use.col]] == values[13] |
                   dat[[use.col]] == values[14] |
                   dat[[use.col]] == values[15] |
                   dat[[use.col]] == values[16] |
                   dat[[use.col]] == values[17] |
                   dat[[use.col]] == values[18] |
                   dat[[use.col]] == values[19] |
                   dat[[use.col]] == values[20] |
                   dat[[use.col]] == values[21] |
                   dat[[use.col]] == values[22] |
                   dat[[use.col]] == values[23] |
                   dat[[use.col]] == values[24] |
                   dat[[use.col]] == values[25] |
                   dat[[use.col]] == values[26] |
                   dat[[use.col]] == values[27] |
                   dat[[use.col]] == values[28] |
                   dat[[use.col]] == values[29] |
                   dat[[use.col]] == values[30] |
                   dat[[use.col]] == values[31], ]
    }
    
    if (n.values == 32) {
      res <- dat[dat[[use.col]] == values[1] |
                   dat[[use.col]] == values[2] |
                   dat[[use.col]] == values[3] |
                   dat[[use.col]] == values[4] |
                   dat[[use.col]] == values[5] |
                   dat[[use.col]] == values[6] |
                   dat[[use.col]] == values[7] |
                   dat[[use.col]] == values[8] |
                   dat[[use.col]] == values[9] |
                   dat[[use.col]] == values[10] |
                   dat[[use.col]] == values[11] |
                   dat[[use.col]] == values[12] |
                   dat[[use.col]] == values[13] |
                   dat[[use.col]] == values[14] |
                   dat[[use.col]] == values[15] |
                   dat[[use.col]] == values[16] |
                   dat[[use.col]] == values[17] |
                   dat[[use.col]] == values[18] |
                   dat[[use.col]] == values[19] |
                   dat[[use.col]] == values[20] |
                   dat[[use.col]] == values[21] |
                   dat[[use.col]] == values[22] |
                   dat[[use.col]] == values[23] |
                   dat[[use.col]] == values[24] |
                   dat[[use.col]] == values[25] |
                   dat[[use.col]] == values[26] |
                   dat[[use.col]] == values[27] |
                   dat[[use.col]] == values[28] |
                   dat[[use.col]] == values[29] |
                   dat[[use.col]] == values[30] |
                   dat[[use.col]] == values[31] |
                   dat[[use.col]] == values[32], ]
    }
    
    if (n.values == 33) {
      res <- dat[dat[[use.col]] == values[1] |
                   dat[[use.col]] == values[2] |
                   dat[[use.col]] == values[3] |
                   dat[[use.col]] == values[4] |
                   dat[[use.col]] == values[5] |
                   dat[[use.col]] == values[6] |
                   dat[[use.col]] == values[7] |
                   dat[[use.col]] == values[8] |
                   dat[[use.col]] == values[9] |
                   dat[[use.col]] == values[10] |
                   dat[[use.col]] == values[11] |
                   dat[[use.col]] == values[12] |
                   dat[[use.col]] == values[13] |
                   dat[[use.col]] == values[14] |
                   dat[[use.col]] == values[15] |
                   dat[[use.col]] == values[16] |
                   dat[[use.col]] == values[17] |
                   dat[[use.col]] == values[18] |
                   dat[[use.col]] == values[19] |
                   dat[[use.col]] == values[20] |
                   dat[[use.col]] == values[21] |
                   dat[[use.col]] == values[22] |
                   dat[[use.col]] == values[23] |
                   dat[[use.col]] == values[24] |
                   dat[[use.col]] == values[25] |
                   dat[[use.col]] == values[26] |
                   dat[[use.col]] == values[27] |
                   dat[[use.col]] == values[28] |
                   dat[[use.col]] == values[29] |
                   dat[[use.col]] == values[30] |
                   dat[[use.col]] == values[31] |
                   dat[[use.col]] == values[32] |
                   dat[[use.col]] == values[33], ]
    }
    
    if (n.values == 34) {
      res <- dat[dat[[use.col]] == values[1] |
                   dat[[use.col]] == values[2] |
                   dat[[use.col]] == values[3] |
                   dat[[use.col]] == values[4] |
                   dat[[use.col]] == values[5] |
                   dat[[use.col]] == values[6] |
                   dat[[use.col]] == values[7] |
                   dat[[use.col]] == values[8] |
                   dat[[use.col]] == values[9] |
                   dat[[use.col]] == values[10] |
                   dat[[use.col]] == values[11] |
                   dat[[use.col]] == values[12] |
                   dat[[use.col]] == values[13] |
                   dat[[use.col]] == values[14] |
                   dat[[use.col]] == values[15] |
                   dat[[use.col]] == values[16] |
                   dat[[use.col]] == values[17] |
                   dat[[use.col]] == values[18] |
                   dat[[use.col]] == values[19] |
                   dat[[use.col]] == values[20] |
                   dat[[use.col]] == values[21] |
                   dat[[use.col]] == values[22] |
                   dat[[use.col]] == values[23] |
                   dat[[use.col]] == values[24] |
                   dat[[use.col]] == values[25] |
                   dat[[use.col]] == values[26] |
                   dat[[use.col]] == values[27] |
                   dat[[use.col]] == values[28] |
                   dat[[use.col]] == values[29] |
                   dat[[use.col]] == values[30] |
                   dat[[use.col]] == values[31] |
                   dat[[use.col]] == values[32] |
                   dat[[use.col]] == values[33] |
                   dat[[use.col]] == values[34], ]
    }
    
    if (n.values == 35) {
      res <- dat[dat[[use.col]] == values[1] |
                   dat[[use.col]] == values[2] |
                   dat[[use.col]] == values[3] |
                   dat[[use.col]] == values[4] |
                   dat[[use.col]] == values[5] |
                   dat[[use.col]] == values[6] |
                   dat[[use.col]] == values[7] |
                   dat[[use.col]] == values[8] |
                   dat[[use.col]] == values[9] |
                   dat[[use.col]] == values[10] |
                   dat[[use.col]] == values[11] |
                   dat[[use.col]] == values[12] |
                   dat[[use.col]] == values[13] |
                   dat[[use.col]] == values[14] |
                   dat[[use.col]] == values[15] |
                   dat[[use.col]] == values[16] |
                   dat[[use.col]] == values[17] |
                   dat[[use.col]] == values[18] |
                   dat[[use.col]] == values[19] |
                   dat[[use.col]] == values[20] |
                   dat[[use.col]] == values[21] |
                   dat[[use.col]] == values[22] |
                   dat[[use.col]] == values[23] |
                   dat[[use.col]] == values[24] |
                   dat[[use.col]] == values[25] |
                   dat[[use.col]] == values[26] |
                   dat[[use.col]] == values[27] |
                   dat[[use.col]] == values[28] |
                   dat[[use.col]] == values[29] |
                   dat[[use.col]] == values[30] |
                   dat[[use.col]] == values[31] |
                   dat[[use.col]] == values[32] |
                   dat[[use.col]] == values[33] |
                   dat[[use.col]] == values[34] |
                   dat[[use.col]] == values[35], ]
    }
  
    gc()
  
    ### Wrap up
  
    return(res)
}
