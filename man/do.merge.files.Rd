% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/do.merge.files.R
\name{do.merge.files}
\alias{do.merge.files}
\title{do.merge.files - Function to merge a list of data.tables (one data.table per 'sample') into a single large data.table.}
\usage{
do.merge.files(dat, remove.duplicates)

do.merge.files(dat, remove.duplicates = TRUE)
}
\arguments{
\item{dat}{NO DEFAULT. List of data.tables (or data.frames)}

\item{remove.duplicates}{DEFAULT = TRUE. Do you want to remove duplicates?}
}
\value{
Returns a combined data.table.
}
\description{
do.merge.files - Function to merge a list of data.tables (one data.table per 'sample') into a single large data.table.
}
\examples{
data.list <- list(
"01_Mock_01" = Spectre::demo.clustered[Sample == '01_Mock_01'],
"02_Mock_02" = Spectre::demo.clustered[Sample == '02_Mock_02']
)
cell.dat <- do.merge.files(dat = data.list, remove.duplicates = TRUE)
head(cell.dat)

}
\references{
\url{https://sydneycytometry.org.au/spectre}.
}
\author{
Thomas M Ashhurst, \email{thomas.ashhurst@sydney.edu.au}
}
