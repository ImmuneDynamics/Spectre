#' run.harmony - dun Harmony alignment on a data.table
#'
#' This function allows you to run the 'Harmony' data alignment algorithm on single cell or cytometry data stored in a data.table
#' 
#' @usage run.harmony()
#' 
#' @param dat NO DEFAULT. A data.table with all of the data you wish to align
#' @param align.cols NO default. The columns you wish to align. For cytometry data, this can be the markers themselves or principle components. For single-cell seq data, principle components are recommended.
#' @param batch.col NO default. The column that denotes the batch or dataset that each cell belongs to
#' @param append.name DEFAULT = '_harmony'. Text that will be appended to the new columns containing aligned data
#' 
#' @return Returns a data.table with aligned data added in new columns.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @examples
#' cell.dat <- run.harmony()
#'
#' @import data.table
#'
#' @export

run.harmony <- function(dat,
                        align.cols,
                        batch.cols,
                        append.name = '_harmony'){
  
  ### Packages
  
      if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
      if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')
      if(!is.element('harmony', installed.packages()[,1])) stop('harmony is required but not installed. You can install harmony by running devtools::install_github("immunogenomics/harmony")')
      
  ### Require packages
  
      require(Spectre)
      require(data.table)
      require(harmony)
  
  ### Data prep
  
      message("run.harmony - preparing data (1/3)")
  
      start.dat <- dat

      dat <- dat[,align.cols, with = FALSE]
      nms <- names(dat)
      dat <- as.matrix(dat)
      
      meta <- data.table()
      meta$CellID <- c(1:nrow(dat))
      meta$CellID <- as.character(meta$CellID)
      meta <- cbind(meta, start.dat[,batch.col, with = FALSE])
      meta <- as_tibble(meta)
      
  ### Run harmony
      
      message("run.harmony - running harmony (2/3)")
      
      hrm.res <- harmony::HarmonyMatrix(dat, meta, batch.col, do_pca = FALSE, verbose=FALSE)
      hrm.res <- as.data.table(hrm.res)
      names(hrm.res) <- paste0(names(hrm.res), append.name)
      hrm.res
      
  ### Final preparation and return
      
      final.res <- cbind(start.dat, hrm.res)
      message("run.harmony - harmony complete, returning data (3/3)")
      
      return(final.res)
}
