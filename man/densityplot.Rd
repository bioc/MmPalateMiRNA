\name{densityplot}
\docType{methods}
\alias{densityplot}
\alias{densityplot-methods}
\alias{densityplot,RGList,missing-method}
\alias{densityplot,list,missing-method}

\title{Density plots of log2 intensity values}

\description{
  Plots the estimated density of log2 intensity values for two-color microarrays
}

\usage{

\S4method{densityplot}{RGList,missing}(
    x,
    channel=c("G", "R"),
    group=NULL,
    subset=NULL, 
    \dots)

\S4method{densityplot}{list,missing}(
    x,
    channel=c("G", "R"),
    group=NULL,
    subset=NULL, 
    \dots)
}

\arguments{
  \item{x}{Either an \code{\link[limma:RGList-class]{RGList}} object,
    or a list containing \code{\link[limma:MAList-class]{MAList}}
    and/or \code{\link[Biobase:NChannelSet-class]{NChannelSet}} objects}
  \item{channel}{The channel to use for calculating distances, one of
    either "G" (green or control channel) or "R" (red or experimental
    channel)}
  \item{group}{An optional character string specifying the name of a
    factor to create separate panel displays, which must be in
    \code{x$genes} (for \code{\link[limma:RGList-class]{RGList}}
    objects)}
  \item{subset}{An optional character vector specifying the which levels
    of \code{group} to use in creating separate panel displays}
  \item{\dots}{arguments to pass to \code{\link[lattice:densityplot]{densityplot}}}
}

\section{Methods}{
  \describe{
    
    \item{\code{signature(x = "RGList", data = "missing")}}{
      For \code{\link[limma:RGList-class]{RGList}} objects, separate panel
      displays 
      can be produced for different types of probes, as determined by the
      \code{group} argument.  
    }
    
    \item{\code{signature(x = "list", data = "missing")}}{
      The method for \code{list} objects is intended
      to work with lists of normalized data sets, as either 
      \code{\link[limma:MAList-class]{MAList}} or 
      \code{\link[Biobase:NChannelSet-class]{NChannelSet}} objects.  This
      method will produce separate panel displays for each normalized data
      set, additionally subsetted by the \code{group} argument if supplied.
      The \code{useOuterStrips} function in the \pkg{latticeExtra} package
    can be used for `outer' strip labels in the latter case.
  }
}}



\references{
  D. Sarkar, R. Parkin, S. Wyman, A. Bendoraite, C. Sather, J. Delrow, A. K. Godwin,
  C. Drescher, W. Huber, R. Gentleman, and M. Tewari.
  Quality assessment and data analysis for microRNA expression arrays.
  Nucleic Acids Res, 37(2):e17, 2009.
}

\seealso{
  \code{\link[MmPalateMiRNA:levelplot]{levelplot}} for pairwise distance
  plots between arrays, 
  \code{\link[MmPalateMiRNA:MADvsMedianPlot]{MADvsMedianPlot}} for median
  absolute deviation versus median plots, and 
  \code{\link[MmPalateMiRNA:MAplot]{MAplot}} for MA plots 
}

\examples{
data(PalateData)
res <- densityplot(PalateData, channel="G", group="probe.type", 
                   subset = c("Other miRNAs",   "MMU miRNAs", "Control"), 
                   col=rep(1:3, each=3), lty=rep(1:3, 3), 
                   key = list(lines=list(col=rep(1:3, each=3), lty=rep(1:3, 3)), 
                     columns=3))
print(res)
}

\keyword{methods}
\keyword{hplot} %% High-Level Plots



