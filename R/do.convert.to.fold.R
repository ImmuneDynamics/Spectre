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
      # x <- read.csv("/Users/thomasa/Google Drive File Stream/My Drive/_Sydney Cytometry/Libraries (synced)/GitHub/Public github/Spectre/scripts/SumTables_Graphing/Output-sumtables/SumTable-CellCounts.csv")
      # as.matrix(names(x))
      # sample.col <- "Sample"
      # group.col <- "Group"
      # ctrl.grp <- "Mock"
      # convert.cols <- names(x)[c(6:20)]
      # log2 = TRUE

  ###

      #x <- dat
      x <- as.data.table(x)

      annot <- x[,!convert.cols, with = FALSE]
      temp <- x[,convert.cols, with = FALSE]

      # annot <- x[-c(convert.cols)]
      # temp <- x[c(convert.cols)]

      ctrl.dat <- x[x[[group.col]] == ctrl.grp,]

      #ctrl.grp <- subset(x, x[group.col] == ctrl.grp)
          #ctrl.grp <- ctrl.grp[,unlist(lapply(ctrl.grp, is.numeric))]

      ctrl.grp.means <- colMeans(ctrl.dat[,convert.cols, with = FALSE])
      as.matrix(ctrl.grp.means)

      fold.raw <- t(t(temp) / ctrl.grp.means)
      fold.raw

      if(log2 == TRUE){
        fold <- log(x = fold.raw, 2)
      }

      if(log2 == FALSE){
        fold <- fold.raw
      }

      fold <- as.data.table(fold)

      res <- cbind(annot, fold)

      # x[,convert.cols, with = FALSE] <- fold

      # x[c(convert.cols)] <- fold
      # x

      return(res)
}
