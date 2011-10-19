\name{checkMVs}
%% \Rdversion{1.1}
\docType{methods}
\alias{checkMVs}
\alias{checkMVs-methods}
\alias{checkMVs,RGList-method}

\title{
  Check an \code{\link[limma:RGList-class]{RGList}} object for missing values
}
\description{
  Checks each of the red and green foreground and background channels in
  an \code{\link[limma:RGList-class]{RGList}} object for missing values.  
}
\usage{
checkMVs(obj)
\S4method{checkMVs}{RGList}(obj)
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


\value{
  Returns a list with the following components
  \item{R.na}{index of missing values in the red channel (\code{obj$R})}
  \item{Rb.na}{index of missing values in the red background channel (\code{obj$Rb})}
  \item{G.na}{index of missing values in the green channel (\code{obj$G})}
  \item{Gb.na}{index of missing values in the green background channel (\code{obj$Gb})}

}

\seealso{
  \code{\link{fixMVs}},  \code{\link{checkOutliers}}, \code{\link{fixOutliers}}
}

\examples{
data(PalateData)
mvs <- checkMVs(PalateData)
}

\keyword{methods}
\keyword{manip} %% data manipulation 
