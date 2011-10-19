\name{fixOutliers}
\alias{fixOutliers}
\title{
  'Fix' an \code{RGList} object with outlying values.  
}
\description{
  Imputes outlying values in one of the red foreground, red background,
  green foreground, or green background matrices of an \code{RGList}
  object created from Miltenyi Biotech miRXplore Microarrays.
  Uses the replicate structure of the array to impute the outlying
  values. Implicit assumption is that only
  one of the four replicated values for a probe is an outlying value.  
}
\usage{
fixOutliers(mat, idx, gene.ids)
}
\arguments{
  \item{mat}{
    One of the red foreground (\code{R}), red background (\code{Rb}),
    green foreground (\code{G}), or green background (\code{Gb})
    matrices in an \code{RGList} object, which contains
    outlying values.}
  \item{idx}{
    Index of outlying values, as returned by the
    \code{\link{checkOutliers}} function.  See examples for usage.}
\item{gene.ids}{
      Vector of gene IDs for each probe. See examples for usage.}
}
\details{
  The function is specific to \code{RGList} objects which were created
  from  Miltenyi Biotech miRXplore Microarrays,
  since it depends on the replicated
  structure of that array (probes spotted in quadruplicate) to impute
  the outlying probe values.  
}
\value{
  Returns a matrix with the outlying probe values imputed.
}


\seealso{
  \code{\link{checkOutliers}}, \code{\link{checkMVs}}, \code{\link{fixMVs}}
}
\examples{
data(PalateData)
outliers <- checkOutliers(PalateData)
PalateData$R <- fixOutliers(PalateData$R, outliers$Rout, PalateData$genes$Gene)
}

\keyword{manip}   %% data manipulation 