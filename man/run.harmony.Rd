% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.harmony.R
\name{run.harmony}
\alias{run.harmony}
\title{run.harmony - dun Harmony alignment on a data.table}
\usage{
run.harmony()
}
\arguments{
\item{dat}{NO DEFAULT. A data.table with all of the data you wish to align}

\item{align.cols}{NO default. The columns you wish to align. For cytometry data, this can be the markers themselves or principle components. For single-cell seq data, principle components are recommended.}

\item{batch.col}{NO default. The column that denotes the batch or dataset that each cell belongs to}

\item{append.name}{DEFAULT = '_aligned'. Text that will be appended to the new columns containing aligned data}

\item{do_pca}{DEFAULT = TRUE. Whether to perform PCA on input matrix.}

\item{npcs}{If doing PCA on input matrix, number of PCs to compute.}

\item{theta}{Diversity clustering penalty parameter. Specify for each
variable in vars_use Default theta=2. theta=0 does not encourage any
diversity. Larger values of theta result in more diverse clusters.}

\item{lambda}{Ridge regression penalty parameter. Specify for each variable
in vars_use.
Default lambda=1. Lambda must be strictly positive. Smaller values result
in more aggressive correction.}

\item{sigma}{Width of soft kmeans clusters. Default sigma=0.1. Sigma scales
the distance from a cell to cluster centroids. Larger values of sigma
result in cells assigned to more clusters. Smaller values of sigma make
soft kmeans cluster approach hard clustering.}

\item{nclust}{Number of clusters in model. nclust=1 equivalent to simple
linear regression.}

\item{tau}{Protection against overclustering small datasets with large ones.
tau is the expected number of cells per cluster.}

\item{block.size}{What proportion of cells to update during clustering.
Between 0 to 1, default 0.05. Larger values may be faster but less accurate}

\item{max.iter.harmony}{Maximum number of rounds to run Harmony. One round
of Harmony involves one clustering and one correction step.}

\item{max.iter.cluster}{Maximum number of rounds to run clustering at each
round of Harmony.}

\item{epsilon.cluster}{Convergence tolerance for clustering round of
Harmony. Set to -Inf to never stop early.}

\item{epsilon.harmony}{Convergence tolerance for Harmony. Set to -Inf to
never stop early.}

\item{plot_convergence}{Whether to print the convergence plot of the
clustering objective function. TRUE to plot, FALSE to suppress. This can be
useful for debugging.}

\item{return_object}{(Advanced Usage) Whether to return the Harmony object
or only the corrected PCA embeddings.}

\item{verbose}{DEFAULT = FALSE. Whether to print progress messages. TRUE to print,
FALSE to suppress.}

\item{reference_values}{(Advanced Usage) Defines reference dataset(s).
Cells that have batch variables values matching reference_values will not
be moved.}

\item{cluster_prior}{(Advanced Usage) Provides user defined clusters for
cluster initialization. If the number of provided clusters C is less than K,
Harmony will initialize K-C clusters with kmeans. C cannot exceed K.}
}
\value{
Returns a data.table with aligned data added in new columns.
}
\description{
This function allows you to run the 'Harmony' data alignment algorithm on single cell or cytometry data stored in a data.table
}
\author{
Thomas M Ashhurst, \email{thomas.ashhurst@sydney.edu.au}
}
