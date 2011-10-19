\name{MAplot}
\docType{methods}
\alias{MAplot}
\alias{MAplot-methods}
\alias{MAplot,MAList-method}
\alias{MAplot,NChannelSet-method}


\title{MA plot}

\description{
Plots of the log2 expression ratios (M values) versus the mean
log2 expression values (A values) for each probe for each array.  
}

\usage{
MAplot(x, \dots)

\S4method{MAplot}{MAList}(
    x,
    \dots)

\S4method{MAplot}{NChannelSet}(
    x,
    \dots)
}

\arguments{
  \item{x}{Either an \code{\link[limma:MAList-class]{MAList}} object 
    or an \code{\link[Biobase:NChannelSet-class]{NChannelSet}} object}
  \item{\dots}{arguments to pass to \code{\link[lattice:xyplot]{xyplot}}}  
}

\section{Methods}{
  \describe{
    \item{\code{signature(x = "MAList")}}{
      \code{M} and \code{A} values are stored as matrices in \code{x}
    }
    \item{\code{signature(x = "NChannelSet")}}{
      \code{M} and \code{A} values are calculated from the \code{R} and
      \code{G} matrices returned by \code{assayData(x)}
    }
}}


\details{
The  so-called "MA" plot can be used to evaluate whether their is a bias
associated with overall intenstity level for each array.  Loess smoothed
regression lines are superimposed on each plot to demonstrate the trend.
}


\seealso{
\code{\link[MmPalateMiRNA:densityplot]{densityplot}} for density plots
of log2 intensity values,
  \code{\link[MmPalateMiRNA:levelplot]{levelplot}} for pairwise distance
  plots between arrays, 
and \code{\link[MmPalateMiRNA:MADvsMedianPlot]{MADvsMedianPlot}} for median
absolute deviation versus median plots.
}

%% examples on reduced (filtered) data ... 
\examples{
data(PalateData)
reducedSet <- filterArray(PalateData, keep=c("MIR", "LET", "POSCON", "CALIB"),
                          frac=1.1, number=3, reps=4)
ndata.quantile <- normalizeBetweenArrays(reducedSet, method="quantile")
res <- MAplot(ndata.quantile)
print(res)
}

\keyword{methods}
\keyword{hplot} %% High-Level Plots


  
