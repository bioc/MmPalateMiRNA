\name{fixMVs}
\alias{fixMVs}

\title{
  'Fix' an \code{RGList} object with missing values.  
}
\description{
  Imputes missing values in one of the red foreground, red background,
  green foreground, or green background matrices of an \code{RGList}
  object created from Miltenyi Biotech miRXplore Microarrays.
  Uses the replicate structure of the array to impute the missing
  values. Implicit assumption is that only
  one of the four replicated values for a probe is an missing value.
}
\usage{
fixMVs(mat, idx, gene.ids)
}

\arguments{
  \item{mat}{
    One of the red foreground (\code{R}), red background (\code{Rb}),
    green foreground (\code{G}), or green background (\code{Gb})
    matrices in an \code{RGList} object, which contains
    missing values.}
  \item{idx}{
    Index of missing values, as returned by the
    \code{\link{checkMVs}} function.  See examples for usage.}
  \item{gene.ids}{
    Vector of gene IDs for each probe. See examples for usage.}
}

\details{
  The function is specific to \code{RGList} objects which were created
  from  Miltenyi Biotech miRXplore Microarrays,
  since it depends on the replicated
  structure of that array (probes spotted in quadruplicate) to impute
  the missing probe values.  
}

\value{
  Returns a matrix with the missing probe values imputed.
}


\seealso{
  \code{\link{checkMVs}}, \code{\link{checkOutliers}},
  \code{\link{fixOutliers}}
  }

\examples{
data(PalateData)
mvs <- checkMVs(PalateData)
PalateData$Rb <- fixMVs(PalateData$Rb, mvs$Rb.na, PalateData$genes$Gene)
}

\keyword{manip}   %% data manipulation 
