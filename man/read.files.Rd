% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.files.R
\name{read.files}
\alias{read.files}
\title{read.files - Function to read data from CSV or FCS files into a list.}
\usage{
read.files(file.loc, file.type, do.embed.file.names, header)
}
\arguments{
\item{file.loc}{DEFAULT = getwd(). What is the location of your files?}

\item{file.type}{DEFAULT = ".csv". What type of files do you want to read. Can be ".csv" or ".fcs".}

\item{files}{DEFAULT = NULL. A vector of selected file names to import.}

\item{nrows}{DEFAULT = NULL. Can specify a numerical target for the number of cells (rows) to be read from each file. Please note, order is random in FCS files.}

\item{do.embed.file.names}{DEFAULT = TRUE. Do you want to embed each row (cell) of each file with the name name?}

\item{header}{DEFAULT = TRUE. Does the first line of data contain column names?}
}
\value{
Returns a list of data.tables -- one per CSV file.
}
\description{
This function allows you to read in sample files (.csv or .fcs) into a list, where each file is saved as a data.table.
}
\examples{
# download sample data
download.file(url='https://github.com/ImmuneDynamics/data/blob/main/msCNS.zip?raw=TRUE', destfile = 'msCNS.zip', mode = 'wb')
unzip(zipfile = 'msCNS.zip')
setwd("msCNS/data")
data.list <- read.files(file.type = ".csv", do.embed.file.names = TRUE)

# return to previous working directory
setwd("../../")

}
\references{
Ashhurst, T. M., et al. (2019). \url{https://www.ncbi.nlm.nih.gov/pubmed/31077106}
}
\author{
Thomas M Ashhurst, \email{thomas.ashhurst@sydney.edu.au}
Felix Marsh-Wakefield, \email{felix.marsh-wakefield@sydney.edu.au}
}
