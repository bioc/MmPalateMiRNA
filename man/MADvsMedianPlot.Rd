\name{MADvsMedianPlot}
\docType{methods}
\alias{MADvsMedianPlot}
\alias{MADvsMedianPlot-methods}
\alias{MADvsMedianPlot,list-method}

\title{Spread vs location of probe intensities}

\description{
Plots of the spread (median absolute deviation)
versus the location (median) of probe intensity levels.
}

\usage{

MADvsMedianPlot(x, \dots)

\S4method{MADvsMedianPlot}{list}(
    x,
    channel=c("G", "R"),
    group=NULL,
    subset=NULL, 
    \dots)

}


\arguments{
  \item{x}{A list containing \code{\link[limma:MAList-class]{MAList}}
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
    \item{\code{signature(x = "list")}}{
      The method for \code{list} objects is intended
      to work with lists of normalized data sets, as either 
      \code{\link[limma:MAList-class]{MAList}} or 
      \code{\link[Biobase:NChannelSet-class]{NChannelSet}} objects.  This
      method will produce separate panel displays for each normalized data
      set, additionally color-coded by the \code{group} argument if supplied.
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
  \code{\link[MmPalateMiRNA:densityplot]{densityplot}} for density plots
  of log2 intensity values, and
  \code{\link[MmPalateMiRNA:MAplot]{MAplot}} for MA plots. 
}

\examples{
data(PalateData)
reducedSet <- filterArray(PalateData, keep=c("MIR", "LET", "POSCON", "CALIB"),
                          frac=1.1, number=3, reps=4)
ndata.none <- normalizeWithinArrays(reducedSet, method="none")
ndata.median <- normalizeWithinArrays(reducedSet, method="median")
ndata.loess <- normalizeWithinArrays(reducedSet, method="loess")
ndata.quantile <- normalizeBetweenArrays(reducedSet, method="quantile")
ndata.all <- list(ndata.none, ndata.median, ndata.loess, 
                  ndata.quantile)
res <- MADvsMedianPlot(ndata.all, channel="R", group="probe.type",                  
                 subset=c("MMU miRNAs", "Other miRNAs", "Control"))
print(res)

}

\keyword{methods}
\keyword{hplot} %% High-Level Plots



