% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.ruv.R
\name{run.ruv}
\alias{run.ruv}
\title{run.ruv - Run the the CytofRUV alignment method.}
\usage{
run.align()
}
\arguments{
\item{dat}{NO DEFAULT.A data.table consisting of the data (both reference and taret) you wish to align}

\item{sample.col}{NO DEFAULT. Column name of the data.table that contains sample names}

\item{batch.col}{NO DEFAULT. Column name of the data.table that contains batch labels}

\item{align.cols}{NO DEFAULT. Vector of column names to align.}

\item{ref.samples}{NO DEFAULT. Character vector of sample names that are used as reference data.}

\item{cluster.col}{DEFAULT = NULL. Column name of the data.table that contains clusters. If NULL, will align data without clusters.}

\item{k}{DEFAULT = K. K value for alignment}

\item{dir}{DEFAULT = getwd(). Directory where the work is being done}

\item{append.name}{DEFAULT = '_ruv'. Text to be appended to end of each new column containing aligned data}
}
\value{
Returns a data.table with aligned data added in new columns.
}
\description{
This function allows you to run the CytofRUV alignment method.
}
\references{
Ashhurst, T. M., et al. (2019). \url{https://www.ncbi.nlm.nih.gov/pubmed/31077106}
}
\author{
Thomas M Ashhurst, \email{thomas.ashhurst@sydney.edu.au}
}
