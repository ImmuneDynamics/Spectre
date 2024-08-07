% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.flowsom.R
\name{run.flowsom}
\alias{run.flowsom}
\title{Run the FlowSOM algorithm}
\usage{
run.flowsom(dat, use.cols, meta.k, xdim, ydim, clust.seed, 
meta.seed, clust.name, meta.clust.name)
}
\arguments{
\item{dat}{NO DEFAULT. data.frame. Input sample.}

\item{use.cols}{NO DEFAULT. Vector of column names to use for clustering.}

\item{xdim}{DEFAULT = 14. Numeric. Number of first level clusters across the x-axis. xdim x ydim = total number of first level FlowSOM clusters.}

\item{ydim}{DEFAULT = 14. Numeric. Number of first level clusters across the y-axis. xdim x ydim = total number of first level FlowSOM clusters.}

\item{meta.k}{DEFAULT = 'auto'. If set to 'auto', then number of metaclusters will be determined automatically. Alternatively, can specify the desired number of metaclusters to create. If set to zero (0), no metaclusters will be created.}

\item{max.meta}{DEFAULT = 20. Only used if meta.k is set to 'auto'. This parameter indicates the maximum number of metaclusters FlowSOM will try out when determining the optimal number of metaclusters for the dataset.}

\item{clust.seed}{DEFAULT = 42 Numeric. Clustering seed for reproducibility.}

\item{meta.seed}{DEFAULT = 42 Numeric. Metaclustering seed for reproducibility.}

\item{clust.name}{DEFAULT = "FlowSOM_cluster". Character. Name of the resulting 'cluster' parameter.}

\item{meta.clust.name}{DEFAULT = "FlowSOM_metacluster". Character. Name of the resulting 'metacluster' parameter.}

\item{mem.ctrl}{DEFAULT = TRUE. Runs gc() (garbage collection) after a number of steps to free up memory that hasn't been released quickly enough.}
}
\description{
Method to run the FlowSOM clustering algorithm.
This function runs FlowSOM on a data.table with cells (rows) vs markers (columns) with new columns for FlowSOM clusters and metaclusters.
Output data will be "flowsom.res.original" (for clusters) and "flowsom.res.meta" (for metaclusters).
Uses the R packages "FlowSOM" for clustering, "flowCore" for handling .fcs files, "Biobase" for creating a flow frame, "data.table" for handling data.table format.
}
\examples{
# Run FlowSOM on demonstration dataset
res <- Spectre::run.flowsom(Spectre::demo.clustered,
                            use.cols = c("NK11_asinh", "CD3_asinh", 
                            "CD45_asinh", "Ly6G_asinh", "CD11b_asinh", 
                            "B220_asinh", "CD8a_asinh", "Ly6C_asinh", 
                            "CD4_asinh"))

}
\author{
Thomas Ashhurst, \email{thomas.ashhurst@sydney.edu.au}
Felix Marsh-Wakefield, \email{felix.marsh-wakefield@sydney.edu.au}
}
