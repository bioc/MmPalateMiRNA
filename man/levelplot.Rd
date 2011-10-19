\name{levelplot}
\docType{methods}
\alias{levelplot}
\alias{levelplot-methods}
%% NOTE - need to include alias for generic???
%% is ok, can still get help pages for both functions (R will prompt for
%% which one ... 
\alias{levelplot,RGList,missing-method}
\alias{levelplot,list,missing-method}

\title{Pairwise distance between arrays}

\description{
  Calculates and plots the pairwise distance between arrays, as measured
  by the median of the absolute differences in log2 intensity values. 
}


\usage{

\S4method{levelplot}{RGList,missing}(
    x,
    channel=c("G", "R"),
    group=NULL,
    subset=NULL, 
    \dots)

\S4method{levelplot}{list,missing}(
    x,
    channel=c("G", "R"),
    order=NULL, 
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
  \item{order}{An optional numeric vector specifying the order of the
    arrays to use in producing the distance plots, i.e. for grouping certain
    arrays together}
  \item{\dots}{arguments to pass to \code{\link[lattice:levelplot]{levelplot}}}
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
      set. 
    }
}}



\references{
  D. Sarkar, R. Parkin, S. Wyman, A. Bendoraite, C. Sather, J. Delrow, A. K. Godwin,
  C. Drescher, W. Huber, R. Gentleman, and M. Tewari.
  Quality assessment and data analysis for microRNA expression arrays.
  Nucleic Acids Res, 37(2):e17, 2009.
}

\seealso{
\code{\link[MmPalateMiRNA:densityplot]{densityplot}} for density plots
of log2 intensity values,
\code{\link[MmPalateMiRNA:MADvsMedianPlot]{MADvsMedianPlot}} for median
absolute deviation versus median plots, and 
\code{\link[MmPalateMiRNA:MAplot]{MAplot}} for MA plots 
}

\examples{
data(PalateData)
res <- levelplot(PalateData[, c(1,5,9,2:4,6:8)], 
                 channel="G", group="probe.type", 
                 subset=c("MMU miRNAs", "Other miRNAs", "Control", "Empty"), 
                 scales = list(rot=c(45, 45))) 
print(res)
}

\keyword{methods}
\keyword{hplot} %% High-Level Plots



