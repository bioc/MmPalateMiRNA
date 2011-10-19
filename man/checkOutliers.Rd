\name{checkOutliers}
\docType{methods}
\alias{checkOutliers}
\alias{checkOutliers-methods}
\alias{checkOutliers,RGList-method}

\title{
  Check \code{\link[limma:RGList-class]{RGList}} object for outlying values
}

\description{
  Checks each of the red and green foreground and background channels in
  an \code{\link[limma:RGList-class]{RGList}} for outlying values.
}

\usage{
checkOutliers(obj)
\S4method{checkOutliers}{RGList}(obj)
}

\arguments{
  \item{obj}{
    An \code{\link[limma:RGList-class]{RGList}} object
}
}

\section{Methods}{
  \describe{
    \item{\code{signature(obj = "RGList")}}{}
}}


\details{
  Detects outliers outside range of mean +/- 2.665 standard deviations.
  Returns the indexes of outlying observations in each channel (R,Rb and G,Gb).
}

\value{
Returns a list with the following components
  \item{Rout}{index of outliers in the red channel (\code{obj$R})}
  \item{Rbout}{index of outliers in the red background channel (\code{obj$Rb})}
  \item{Gout}{index of outliers in the green channel (\code{obj$G})}
  \item{Gbout}{index of outliers in the green background channel (\code{obj$Gb})}
}


\seealso{
  \code{\link{fixOutliers}}, \code{\link{checkMVs}}, \code{\link{fixMVs}}
}

\examples{
data(PalateData)
outliers <- checkOutliers(PalateData)
}

\keyword{methods}
\keyword{manip}   %% data manipulation 


