#' write.sumtables - Export generated data as a .csv file
#'
#' This function summarises generated results and exports them as a .csv file.
#' It can include the calculation of median fluorescence intensity (or equivalent), proportion (%) or cell counts for clusters/subsets of interest.
#' Makes use of the packages 'data.table' and 'tidyr' to handle the data.
#'
#' @param dat NO DEFAULT. A data.table containing cells (rows) vs features/markers (columns). One column must represent sample names, and another must represent populations/clusters.
#' @param sample.col NO DEFAULT. Character. Name of the sample column (e.g. "Sample").
#' @param pop.col NO DEFAULT. Character. Name of the population/cluster column (e.g. "Population", "Cluster").
#' @param measure.col NO DEFAULT. A character or numeric vector indicating the columns to be measured (e.g. cellular columns -- c("CD45", "CD3e") etc).
#'
#' @param annot.col DEFAULT = NULL. A character or numeric vector indicating the columns to be included as annotation columns (e.g. cellular columns -- c("Batch", "Group") etc). If groups are present, this column must be present here also.
#' @param group.col DEFAULT = NULL. Character. Which column represent groups (e.g. "Group"). This is for dividing data within the function. If you wish for groups to be included in the annotations, you will need to also include it in the annot.col argument.

#' @param do.proportions DEFAULT = TRUE. Do you wish to create cell proportion and/or cell count results?
#' @param cell.counts DEFAULT = NULL. If you wish to generate cell.count results, a vector of cell counts (e.g. c(1000, 1500, 2439,)) representing the cell counts in each of the samples. Must be entered in the order the that unique sample names appear in the dataset.
#' @param do.mfi.per.sample DEFAULT = TRUE. Do you wish to generate MFI data (markers vs clusters) for each sample?
#' @param do.mfi.per.marker DEFAULT = FALSE. Do you wish to generate MFI data (sample vs clusters) for each marker?
#'
#' @param perc.pos.markers DEFAULT = NULL. A vector of column names of calculating percent positive summary stats.
#' @param perc.pos.cutoff DEFAULT = NULL. A vector of 'positive' cut-off values for the markers defined in perc.pos.markers. Must be in same order.
#'
#' @param mfi.type DEFAULT = "median". Can be "median" or "mean". Defines the type of function for calculating MFI data.
#' @param path DEFAULT = getwd(). Defines the directory for write CSV files.
#' 
#' @usage make.sumtables(dat, sample.col, pop.col, measure.col, annot.col, group.col, do.proportions, cell.counts, do.mfi.per.sample, do.mfi.per.marker, perc.pos.markers, perc.pos.cutoff, mfi.type, path)
#'
#' @examples
#' # Calculate and export results from demonstration data
#' Spectre::write.sumtables(dat = Spectre::demo.clustered,
#'                          sample.col = "Sample",
#'                          pop.col = "FlowSOM_metacluster",
#'                          measure.col = c(2,5:6,8:9,11:13,16:19,21:30,32),
#'                          annot.col = c(33:34,36:37),
#'                          group.col = "Group",
#'                          cell.counts = c(rep(2.0e+07, 6), rep(1.8e+07, 6)),
#'                          do.mfi.per.marker = TRUE,
#'                          perc.pos.markers = c("BV711.SCA.1","APC.BrdU"),
#'                          perc.pos.cutoff = c(580, 450)
#'                          )
#'
#' @author
#' Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#' @export

write.sumtables <- function(dat,
                           sample.col, # name
                           pop.col, # name
                           measure.col, # numbers

                           ## Optional extras -- with implications
                           annot.col = NULL, # numbers # takes the first value from each 'sample' -- as these are annotations that are entered per 'cell', but should be equivalent for each sample
                           group.col = NULL, # name

                           ## Functions
                           do.proportions = TRUE,
                           cell.counts = NULL,
                           do.mfi.per.sample = TRUE,
                           do.mfi.per.marker = FALSE, # coming soon

                           ## Functions for percentage positive calculations
                           perc.pos.markers = NULL,
                           perc.pos.cutoff = NULL,

                           ## Other defaults
                           mfi.type = "median",
                           path = getwd())
{

  message("The 'write.sumtables' function has been depreciated. Please use 'create.sumtable' instead")

}

