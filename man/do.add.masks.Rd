% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/do.add.masks.R
\name{do.add.masks}
\alias{do.add.masks}
\title{do.add.masks}
\usage{
do.add.masks(dat, mask.dir, mask.pattern, mask.label,correct.extent = TRUE, 
flip.y = TRUE, value.modifier = 65535, HDF5 = FALSE, array = FALSE)
}
\arguments{
\item{dat}{NO DEFAULT. A list of spatial data objects}

\item{mask.dir}{NO DEFAULT. Directory of mask files.}

\item{mask.pattern}{NO DEFAULT. A character pattern that identifies the mask type (e.g. '_mask')}

\item{mask.label}{NO DEFAULT. What do you want to call the mask}

\item{correct.extent}{DEFAULT = TRUE}

\item{flip.y}{DEFAULT = TRUE}

\item{value.modifier}{DEFAULT = 65535}

\item{HDF5}{DEFAULT = FALSE. Can read in HDF5 mask files if desired (advanced use).}

\item{array}{DEFAULT = FALSE. If the mask TIFF is an array (advanced use).}
}
\description{
do.add.masks
}
\references{
\url{https://github.com/ImmuneDynamics/Spectre}.
}
\author{
Thomas M Ashhurst, \email{thomas.ashhurst@sydney.edu.au}
}
