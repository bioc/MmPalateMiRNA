\name{PalateData}
\alias{PalateData}
\docType{data}
\title{
  Murine Secondary Palate Development miRNA Expression Data
}
\description{
  This data set contains two-color miRNA microarray expression data
  obtained from mouse embryonic tissue during gestational days (GD) 12,
  13, and 14,  which represents the critical period of palate
  development in the mouse.
}
\usage{data(PalateData)}
\format{
  The data are in the format of an \code{"\linkS4class{RGList}"}, 
  which in this case is a list with the following 9 elements:
  \describe{
    \item{\code{R}}{matrix of dimension 6336 x 9 which contains the red
      channel foreground measurements}
    \item{\code{G}}{matrix of dimension 6336 x 9 which contains the green
      channel foreground measurements}
    \item{\code{Rb}}{matrix of dimension 6336 x 9 which contains the red
      channel background measurements}
    \item{\code{Gb}}{matrix of dimension 6336 x 9 which contains the green
      channel background measurements}
    \item{\code{source}}{source of the images, here "imagene"}
    \item{\code{Field.Dimensions}}{numeric vector giving the field
      dimensions of the array (Metarows, Metacols, Rows and Cols)}
    \item{\code{weights}}{matrix of dimension 6336 x 9 which contains
      the quality weights associated with each spot on the arrays}
    \item{\code{printer}}{list containing information on the process
      used to print the spots on the arrays (number of grid rows /
      columns and number of spot rows / columns per grid - coincides
      with \code{Field.Dimensions})}
    \item{\code{genes}}{ A \code{data.frame} containing information on
      each probe. Has the following columns:
      \describe{
	\item{\code{Field}}{field position for the probe}
	\item{\code{Meta Row}}{meta row position for the probe}
	\item{\code{Meta Column}}{meta column position for the probe}
	\item{\code{Row}}{row position for the probe}
	\item{\code{Column}}{column position for the probe}
	\item{\code{Gene ID}}{unique gene identifier provided by Miltenyi Biotec}
	\item{\code{ID}}{unique probe identifier constructed by
	  concatentating the "Gene ID" with "Meta Row", "Meta Column",
	  "Row", and "Column" information}
	\item{\code{Name}}{name of the microRNA}
	\item{\code{Name.stem}}{base name of the microRNA}
	\item{\code{probe.type}}{type of probe, "MMU miRNAs", "Other
	  miRNAs", "Control", "Empty", and "Other"}
      }
    }
  }
}


\details{
  RNA samples were isolated from mouse embryonic orofacial
  tissues (GD-12 - GD-14) and fluorescently labeled with Hy5 (red).
  Control samples (miRXplore Universal Reference) were labeled with 
  Hy3 (green).  The two sets of samples were hybridized to miRXplore
  Microarrays (Miltenyi Biotec) sing the a-Hyb Hybridization Station
  (Miltenyi Biotec). Probes for a
  total of 1336 mature miRNAs (from human, mouse, rat and virus),
  including positive control and calibration probes, were spotted in
  quadruplicate on each microarray. Each array included probes for 588
  murine miRNAs.
}

\source{

  P. Mukhopadhyay, G. Brock, V. Pihur, C. Webb, M.M. Pisano, and R.M. Greene.
  Developmental microRNA expression profiling of murine embryonic orofacial tissue.
  Birth Defects Res A Clin Mol Teratol, 88(7):511-34, 2010.
  
}

\examples{
data(PalateData)
}
\keyword{datasets}
