\name{clustPlot}
\alias{clustPlot}

\title{
  Plot expression profiles
}
\description{
  Produces plots of clustered expression profiles, with seperate plots
  for each cluster.  The average expression profile for each cluster is
  superimposed as well.
}
\usage{
clustPlot(cl, mat, nrow, ncol)
}

\arguments{
  \item{cl}{integer vector giving the cluster membership for each item }

  \item{mat}{matrix of values to be plotted }

  \item{nrow}{number of rows to use for plotting }

  \item{ncol}{number of columns to use for plotting }

}
\details{
  The figure region will be subdivided into \code{nrow} by \code{ncol}
  separate plots, using \code{mfrow}.  The average expression profile
  and the number of genes belonging to each cluster is superimposed on
  each of the plots.
}

\references{
  G.N. Brock, V. Pihur, S. Datta, and S. Datta.
  clValid, an R package for cluster validation.
  Journal of Statistical Software, 25, 2008.
}


\seealso{
  See the package vignette for illustration on usage
}
\examples{

## generate some fake data and cluster
set.seed(101)
mat <- matrix(rnorm(500), nrow=100, ncol=5)
clusts <- hclust(dist(mat))
cl <- cutree(clusts, 6)
clustPlot(cl, mat, 3, 2)

}

\keyword{ cluster }
