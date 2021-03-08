#' do.add.cols
#'
#' @usage do.add.cols(dat, base.col, add.dat, add.by, rmv.ext)
#'
#' @param dat NO DEFAULT. A data.table (or data.frame) containing the data to have new values added to.
#' @param base.col NO DEFAULT. Name of the column containing values that new values will be matched to.
#' @param add.dat NO DEFAULT. A data table of new values to embed as a new columns, with one column containing values used for matching with the target data.table.
#' @param add.by NO DEFAULT. Character, name of the column in add.dat that is used for matching to the taret dataset.
#' @param rmv.ext DEFAULTS TO FALSE. Logical, can be TRUE or FALSE. Removes a ".csv" or ".fcs" extension from a the 'match.to' vector -- especially useful if 'match.to' is a list of sample names that end in .csv or .fcs.
#' @param mem.ctrl DEFAULT = TRUE. Runs gc() (garbage collection) after a number of steps to free up memory that hasn't been released quickly enough.
#'
#' @return Returns the data.table with new columns embedded.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#'
#' @examples
#' add.dt <- data.table('Files' = unique(Spectre::demo.asinh$FileName),
#'                      'SampleNumber' = c(1:12))
#' 
#' cell.dat <- do.add.cols(dat = Spectre::demo.asinh,
#'                         base.col = "FileName",
#'                         add.dat = add.dt,
#'                         add.by = "Files")
#'
#' @import data.table
#'
#' @export

do.add.cols <- function(dat, # the list of dataframes (samples) where each dataframe will have the columns embedded
                         base.col, # column name in the actual dataset

                         add.dat,
                         add.by,
                         rmv.ext = FALSE,
                         mem.ctrl = TRUE,
                         show.status = TRUE)
{

  ## Packages
      if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
      if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')

  ### Require packages
      require(Spectre)
      require(data.table)

  ### Test data
  
      ## Normal testing
          # dat <- as.data.table(Spectre::demo.start)
          # base.col <- "FileName"
          # 
          # add.dat <- data.table(Filename = unique(dat$FileName),
          #                       NewColA = c('A','B','C','D','E','F','G','H','I','J','K', 'L'),
          #                       NewColB = c(1,2,3,4,5,6,7,8,9,10,11,12))
          # 
          # add.by <- "Filename"
          # rmv.ext = FALSE
          # 
          # mem.ctrl = TRUE
          # 
          # set.seed(42)
          # rw <- sample(nrow(dat))
          # dat <- dat[rw,]

      ## Testing missing values in add.dat
          # add.dat <- data.table(Filename = unique(dat$FileName)[c(1:11)],
          #                       NewColA = c('A','B','C','D','E','F','G','H','I','J','K'),
          #                       NewColB = c(1,2,3,4,5,6,7,8,9,10,11))
   
  ### Checks

      dat <- as.data.table(dat)
      add.dat <- as.data.table(add.dat)

      if(!(base.col %in% colnames(dat))){
        stop(paste0("Your entry '", base.col, "' (base.col) was not found in your dataset (dat). Please make sure you have correctly entered the name of the column."))
      }

      if(!(add.by %in% colnames(add.dat))){
        warning(paste0("Your entry '", add.by, "' (add.by) was not found in your dataset (add.dat). Please make sure you have correctly entered the name of the column."))
      }

      dat.names <- names(dat)
      add.dat.names <- names(add.dat)

      pos <- match(add.by,add.dat.names)
      check.dup <- c(dat.names, add.dat.names[-pos])

      if(any(duplicated(check.dup))){
        stop("You have duplciate column names in your 'dat' and 'add.dat' entries. The 'base.col', and 'add.by' column names may be identical, but all other columns should unique.")
      }

  ### Prep

      if(rmv.ext == TRUE){
        
        if(show.status == TRUE){
          message("Removing '.csv' or '.fcs' extension")
        }
        
          if(is.numeric(add.dat[[add.by]]) == FALSE){
            temp <- add.dat[[add.by]]
            temp <- gsub("*.csv", "", temp)
            temp <- gsub("*.fcs", "", temp)
            add.dat[[add.by]] <- temp
          }
      }

      if(mem.ctrl == TRUE){
        gc()
      }

  ### Mapping data
      if(show.status == TRUE){
        message("Step 1/3. Mapping data")
      }

      added.names <- names(add.dat[,setdiff(names(add.dat),add.by),with = FALSE])
      
      add.dat <- cbind(add.dat[,add.by,with = FALSE],
                       add.dat[,setdiff(names(add.dat),add.by),with = FALSE]
                       )
      
      # p1 <- add.dat[,add.by,with = FALSE]
      # p2 <- add.dat[,setdiff(names(add.dat),add.by),with = FALSE]
      # 
      # added.names <- names(p2)
      # add.dat <- cbind(p1, p2)
      # 
      # rm(p1)
      # rm(p2)

      if(mem.ctrl == TRUE){
        gc()
      }

      names(add.dat)[c(1)] <- base.col

      dat[[base.col]] <- as.factor(dat[[base.col]])
      add.dat[[base.col]] <- as.factor(add.dat[[base.col]])

      if(class(dat[[base.col]]) != class(add.dat[[base.col]])){
        warning(paste0("The column '", base.col, "' in your main dataset", " is ", class(dat[[base.col]]), ", whereas the column '", add.by, "' in the 'to add' data is ", class(add.dat[[base.col]]), ". You may need to adjust their types to ensure they are the same. If you are matching based on character values (filenames, popualtion names etc, then both need to be character."))
      }

  ### Mergind data
      if(show.status == TRUE){
        message("Step 2/3. Merging data")
      }

          ## dat <- merge(dat, add.dat, by = base.col, sort = FALSE, all.x = TRUE)
          dat <- data.table::merge.data.table(dat, add.dat, by = base.col, sort = FALSE, all.x = TRUE)
          
          rm(add.dat)

          if(mem.ctrl == TRUE){
            gc()
          }

          dat <- as.data.table(dat)
          dat <- dat[,c(dat.names,added.names), with = FALSE]

          if(mem.ctrl == TRUE){
            gc()
          }

  ### Returning data
      if(show.status == TRUE){
        message("Step 3/3. Returning data")
      }
      
      return(dat)
}




#' do.embed.columns
#'
#' @usage do.embed.columns(dat, base.col, add.dat, add.by, rmv.ext)
#'
#' @param dat NO DEFAULT. A data.table (or data.frame) containing the data to have new values added to.
#' @param base.col NO DEFAULT. Name of the column containing values that new values will be matched to.
#' @param add.dat NO DEFAULT. A data table of new values to embed as a new columns, with one column containing values used for matching with the target data.table.
#' @param add.by NO DEFAULT. Character, name of the column in add.dat that is used for matching to the taret dataset.
#' @param rmv.ext DEFAULTS TO FALSE Logical, can be TRUE or FALSE. Removes a ".csv" or ".fcs" extension from a the 'match.to' vector -- especially useful if 'match.to' is a list of sample names that end in .csv or .fcs.
#' @param mem.ctrl DEFAULT = TRUE. Runs gc() (garbage collection) after a number of steps to free up memory that hasn't been released quickly enough.
#'
#' @return Returns the data.table with new columns embedded..
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#'
#' @examples
#'
#' @import data.table
#'
#' @export

do.embed.columns <- function(dat, # the list of dataframes (samples) where each dataframe will have the columns embedded
                             base.col, # column name in the actual dataset

                             add.dat,
                             add.by,
                             rmv.ext = FALSE,
                             mem.ctrl = TRUE)
{
  
  message("Warning: Please note, the 'do.embed.columns' function has been depreciated. Please use 'do.add.cols' instead")
  
#   ## Packages
#   if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
#   if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')
# 
#   ### Require packages
#   require(Spectre)
#   require(data.table)
# 
#   ### Test data
#   # dat <- as.data.table(Spectre::demo.start)
#   # base.col <- "FileName"
#   #
#   # add.dat <- data.table(Filename = unique(dat$FileName),
#   #                       NewColA = c('A','B','C','D','E','F','G','H','I','J','K', 'L'),
#   #                       NewColB = c(1,2,3,4,5,6,7,8,9,10,11,12))
#   #
#   # add.dat <- data.table(Filename = unique(dat$FileName)[c(1:11)],
#   #                       NewColA = c('A','B','C','D','E','F','G','H','I','J','K'),
#   #                       NewColB = c(1,2,3,4,5,6,7,8,9,10,11))
#   #
#   # add.by <- "Filename"
#   # rmv.ext = TRUE
# 
#   ### Checks
# 
#   dat <- as.data.table(dat)
#   add.dat <- as.data.table(add.dat)
# 
#   if(!(base.col %in% colnames(dat))){
#     stop(paste0("Your entry '", base.col, "' (base.col) was not found in your dataset (dat). Please make sure you have correctly entered the name of the column."))
#   }
# 
#   if(!(add.by %in% colnames(add.dat))){
#     warning(paste0("Your entry '", add.by, "' (add.by) was not found in your dataset (add.dat). Please make sure you have correctly entered the name of the column."))
#   }
# 
#   dat.names <- names(dat)
#   add.dat.names <- names(add.dat)
# 
#   pos <- match(add.by,add.dat.names)
#   check.dup <- c(dat.names, add.dat.names[-pos])
# 
#   if(any(duplicated(check.dup))){
#     stop("You have duplciate column names in your 'dat' and 'add.dat' entries. The 'base.col', and 'add.by' column names may be identical, but all other columns should unique.")
#   }
# 
#   ### Prep
# 
#   if(rmv.ext == TRUE){
#     message("Removing '.csv' or '.fcs' extension")
# 
#     if(is.numeric(add.dat[[add.by]]) == FALSE){
#       temp <- add.dat[[add.by]]
#       temp <- gsub("*.csv", "", temp)
#       temp <- gsub("*.fcs", "", temp)
#       add.dat[[add.by]] <- temp
#     }
#   }
# 
#   if(mem.ctrl == TRUE){
#     gc()
#   }
# 
#   ### Mapping data
# 
#   dat <- as.data.table(dat)
#   add.dat <- as.data.table(add.dat)
# 
#   message("Step 1/3. Mapping data")
#   p1 <- add.dat[,add.by,with = FALSE]
#   p2 <- add.dat[,setdiff(names(add.dat),add.by),with = FALSE]
# 
#   added.names <- names(p2)
#   add.dat <- cbind(p1, p2)
# 
#   rm(p1)
#   rm(p2)
# 
#   if(mem.ctrl == TRUE){
#     gc()
#   }
# 
#   names(add.dat)[c(1)] <- base.col
# 
#   dat[[base.col]] <- as.factor(dat[[base.col]])
#   add.dat[[base.col]] <- as.factor(add.dat[[base.col]])
# 
#   if(class(dat[[base.col]]) != class(add.dat[[base.col]])){
#     warning(paste0("The column '", base.col, "' in your main dataset", " is ", class(dat[[base.col]]), ", whereas the column '", add.by, "' in the 'to add' data is ", class(add.dat[[base.col]]), ". You may need to adjust their types to ensure they are the same. If you are matching based on character values (filenames, popualtion names etc, then both need to be character."))
#   }
# 
#   ### Mergind data
# 
#   message("Step 2/3. Merging data")
#   res.1 = merge(dat, add.dat, by = base.col, sort = FALSE, all.x = TRUE)
#   res.1 <- as.data.table(res.1)
#   res.2 <- res.1[,c(dat.names,added.names), with = FALSE]
# 
#   if(mem.ctrl == TRUE){
#     gc()
#   }
# 
#   ### Returning data
# 
#   message("Step 3/3. Returning data")
#   return(res.2)
}
