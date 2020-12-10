#' do.filter - filtering data.table using multiple match values
#'
#' This function allows filtering of a data.table using multiple match values -- all cells that contain any of the match values will be filtered, and provided in a new data.table.
#'
#' @param dat NO DEFAULT. A data.table
#' @param use.col DEFAULT = NULL. The column to use for re-ordering
#' @param values DEFAULT = NULL. A vector of values to use for filtering -- all cells that contain any of the match values will be filtered, and provided in a new data.table. Up to 20 match values can be provided
#' 
#' @usage do.filter(dat, use.col, values)
#'
#' @author
#' Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://github.com/ImmuneDynamics/Spectre}.
#'
#' @import data.table
#' 
#' @export

do.filter <- function(dat,
                      use.col,
                      values){
  
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
      
      if(n.values > 20){
        stop('do.filter function cannot take more than 20 target values')
      }
      
  ### Filtering
      
      if(n.values == 1){
        res <- dat[dat[[use.col]] == values[1],]
      }
      
      if(n.values == 2){
        res <- dat[dat[[use.col]] == values[1] |
                     dat[[use.col]] == values[2],
                   ]
      }
      
      if(n.values == 3){
        res <- dat[dat[[use.col]] == values[1] |
                     dat[[use.col]] == values[2] |
                     dat[[use.col]] == values[3],
                   ]
      }
      
      if(n.values == 4){
        res <- dat[dat[[use.col]] == values[1] |
                     dat[[use.col]] == values[2] |
                     dat[[use.col]] == values[3] |
                     dat[[use.col]] == values[4],
                   ]
      }
      
      if(n.values == 5){
        res <- dat[dat[[use.col]] == values[1] |
                     dat[[use.col]] == values[2] |
                     dat[[use.col]] == values[3] |
                     dat[[use.col]] == values[4] |
                     dat[[use.col]] == values[5],
                   ]
      }
      
      if(n.values == 6){
        res <- dat[dat[[use.col]] == values[1] |
                     dat[[use.col]] == values[2] |
                     dat[[use.col]] == values[3] |
                     dat[[use.col]] == values[4] |
                     dat[[use.col]] == values[5] |
                     dat[[use.col]] == values[6],
                   ]
      }
      
      if(n.values == 7){
        res <- dat[dat[[use.col]] == values[1] |
                     dat[[use.col]] == values[2] |
                     dat[[use.col]] == values[3] |
                     dat[[use.col]] == values[4] |
                     dat[[use.col]] == values[5] |
                     dat[[use.col]] == values[6] |
                     dat[[use.col]] == values[7],
                   ]
      }
      
      if(n.values == 8){
        res <- dat[dat[[use.col]] == values[1] |
                     dat[[use.col]] == values[2] |
                     dat[[use.col]] == values[3] |
                     dat[[use.col]] == values[4] |
                     dat[[use.col]] == values[5] |
                     dat[[use.col]] == values[6] |
                     dat[[use.col]] == values[7] |
                     dat[[use.col]] == values[8],
                   ]
      }
      
      if(n.values == 9){
        res <- dat[dat[[use.col]] == values[1] |
                     dat[[use.col]] == values[2] |
                     dat[[use.col]] == values[3] |
                     dat[[use.col]] == values[4] |
                     dat[[use.col]] == values[5] |
                     dat[[use.col]] == values[6] |
                     dat[[use.col]] == values[7] |
                     dat[[use.col]] == values[8] |
                     dat[[use.col]] == values[9],
                   ]
      }
      
      if(n.values == 10){
        res <- dat[dat[[use.col]] == values[1] |
                     dat[[use.col]] == values[2] |
                     dat[[use.col]] == values[3] |
                     dat[[use.col]] == values[4] |
                     dat[[use.col]] == values[5] |
                     dat[[use.col]] == values[6] |
                     dat[[use.col]] == values[7] |
                     dat[[use.col]] == values[8] |
                     dat[[use.col]] == values[9] |
                     dat[[use.col]] == values[10],
                   ]
      }
      
      if(n.values == 11){
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
                     dat[[use.col]] == values[11],
                   ]
      }
      
      if(n.values == 12){
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
                     dat[[use.col]] == values[12],
                   ]
      }
      
      if(n.values == 13){
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
                     dat[[use.col]] == values[13],
                   ]
      }
      
      if(n.values == 14){
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
                     dat[[use.col]] == values[14],
                   ]
      }
      
      if(n.values == 15){
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
                     dat[[use.col]] == values[15],
                   ]
      }
      
      if(n.values == 16){
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
                     dat[[use.col]] == values[16],
                   ]
      }
      
      if(n.values == 17){
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
                     dat[[use.col]] == values[17],
                   ]
      }
      
      if(n.values == 18){
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
                     dat[[use.col]] == values[18],
                   ]
      }
      
      if(n.values == 19){
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
                     dat[[use.col]] == values[19],
                   ]
      }
      
      if(n.values == 20){
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
                     dat[[use.col]] == values[20],
                   ]
      }
  
      gc()
  
  ### Wrap up
  
      return(res)
  
}

