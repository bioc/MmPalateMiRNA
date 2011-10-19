\name{filterArray}
\docType{methods}
\alias{filterArray}
\alias{filterArray-methods}
\alias{filterArray,RGList-method}

\title{
  Filter an \code{\link[limma:RGList-class]{RGList}} object to remove probes
}

\description{
  Filters an \code{\link[limma:RGList-class]{RGList}} object to remove
  probes with foreground 
  intensities not sufficiently above the background intensity.
  Additionally can filter probes based on character strings, to remove
  e.g. control probes.  
}

\usage{
filterArray (obj, \dots)
\S4method{filterArray}{RGList}(obj, keep, frac, number, reps)
}

\arguments{
  \item{obj}{
    An \code{\link[limma:RGList-class]{RGList}} object
  }
  \item{keep}{
    Character vector to be used as a text filter. Only gene names (as
    contained in \code{obj$genes$Name}) which contain these text strings
    will be retained. 
  }
  \item{frac}{
    Fraction to use as a background filter.  Only those probes with
    foreground values (both red and green) greater than 'frac' times the
    background  values pass the filter.
  }
  \item{number}{
    The number of samples required to pass the background filter for
    each probe. 
  }
  \item{reps}{
    The number of replicates for each probe required to pass the
    filtering step.
  }
  \item{\dots}{allows additional arguments to be passed to specific
    methods
  }
}

\section{Methods}{
  \describe{
    \item{\code{signature(obj = "RGList")}}{}
    }}

\value{
  Returns an \code{\link[limma:RGList-class]{RGList}} object identical
  in structure to the input 
  object, but with reduced dimension according to the filtering steps.
}


\seealso{
  \code{\link{checkMVs}},  \code{\link{fixMVs}},
  \code{\link{checkOutliers}}, \code{\link{fixOutliers}}
}

\examples{
data(PalateData)
reducedSet <- filterArray(PalateData, keep=c("MIR", "LET", "POSCON", "CALIB"),
                          frac=1.1, number=3, reps=4) 
}

\keyword{methods}
\keyword{manip} %% data manipulation 