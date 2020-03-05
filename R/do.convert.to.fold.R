#' do.convert.to.fold
#'
#' @export

do.convert.to.fold <- function(x,
                            sample.col,
                            group.col,
                            ctrl.grp,
                            convert.cols,
                            log2 = TRUE)

{

  ### Test data
      # x <- read.csv("/Users/Tom/Google Drive (t.ashhurst@centenary.org.au)/_Sydney Cytometry/2019_Synced/GitHub/Public github/Spectre/worked_example/Output_CAPX/Output-SumTables/Output_CellNums/SumTable_CellsPerTissue.csv")
      # x$X <- NULL
      # as.matrix(names(x))
      # sample.col <- "Sample"
      # group.col <- "Group"
      # ctrl.grp <- "Mock"
      # convert.cols <- c(3:42)
      # log2 = TRUE

  ###
      annot <- x[-c(convert.cols)]
      temp <- x[c(convert.cols)]

      ctrl.grp <- subset(x, x[group.col] == ctrl.grp)
      #ctrl.grp <- ctrl.grp[,unlist(lapply(ctrl.grp, is.numeric))]

      ctrl.grp.means <- colMeans(ctrl.grp[c(convert.cols)])
      as.matrix(ctrl.grp.means)

      fold.raw <- t(t(temp) / ctrl.grp.means)
      fold.raw

      if(log2 == TRUE){
        fold <- log(x = fold.raw, 2)
      }

      if(log2 == FALSE){
        fold <- fold.raw
      }

      x[c(convert.cols)] <- fold
      x

      return(x)
}
