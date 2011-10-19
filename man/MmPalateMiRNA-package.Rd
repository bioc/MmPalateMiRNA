\name{MmPalateMiRNA-package}
\alias{MmPalateMiRNA-package}
\alias{MmPalateMiRNA}
\docType{package}
\title{
  R package compendium for the analysis of murine palate two-color miRNA expression data
}
\description{
  R package compendium for the analysis of two-color miRNA expression
  data, during the period of murine embryonic palate development
  (gestational days (GD) 12, 13, and 14).  Samples were hybridized to
  Miltenyi Biotech miRXplore Microarrays.
  The compendium covers a wide range of
  steps which occur in a  typical miRNA microarray data analysis, including
  pre-processing, normalization, differential expression analysis,
  clustering, target identification, and gene-set enrichment analysis.  
}
\details{
\tabular{ll}{
Package: \tab MmPalateMiRNA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-09-14\cr
License: \tab LGPL\cr
LazyLoad: \tab yes\cr
}
The package contains several functions which are helpful during the
pre-processing steps of array data, which are specific to
\code{\link[limma:RGList-class]{RGList}} objects and, in the case of the
\code{\link{fixOutliers}} and \code{\link{fixMVs}} functions, depend on
the replicated structure of  Miltenyi Biotech miRXplore
Microarrays.  Additionally, methods are available to produce diagnostic plots
for \code{\link[limma:RGList-class]{RGList}} objects and lists of
normalized data sets (\code{\link[limma:MAList-class]{MAList}} and /
or \code{\link[Biobase:NChannelSet-class]{NChannelSet}} objects), which
build on the generic functions in \code{\link[lattice:lattice-package]{lattice}}.
Lastly, the main focus of the package is the package
vignette "MmPalateMiRNA", which contains an extended example covering
the typical steps in an miRNA microarray data analysis.
}
\author{
Guy Brock, Partha Mukhopadhyay, Vasyl Pihur, Bob Green, M. Michele Pisano
Maintainer: Guy Brock <guy.brock@louisville.edu>
}
\references{

  P. Mukhopadhyay, G. Brock, V. Pihur, C. Webb, M.M. Pisano, and R.M. Greene.
  Developmental microRNA expression profiling of murine embryonic orofacial tissue.
  Birth Defects Res A Clin Mol Teratol, 88(7):511-34, 2010.

  R. Gentleman.
  Reproducible research: a bioinformatics case study.
  Stat Appl Genet Mol Biol, 4:Article2, 2005.

  D. Sarkar, R. Parkin, S. Wyman, A. Bendoraite, C. Sather, J. Delrow, A. K. Godwin,
  C. Drescher, W. Huber, R. Gentleman, and M. Tewari.
  Quality assessment and data analysis for microRNA expression arrays.
  Nucleic Acids Res, 37(2):e17, 2009.
  
}

\keyword{ package }
\keyword{ manip }


