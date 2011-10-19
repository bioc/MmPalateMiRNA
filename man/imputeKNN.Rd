\name{imputeKNN}
\alias{imputeKNN}

\title{
  Impute missing values 
}
\description{
  Imputes missing values in a data matrix using the K-nearest neighbor
  algorithm.   
}
\usage{
imputeKNN(data, k = 10, distance = "euclidean", rm.na = TRUE, rm.nan =
TRUE, rm.inf = TRUE )
}

\arguments{
  \item{data}{
    a data matrix
}
  \item{k}{
    number of neighbors to use
}
  \item{distance}{
    distance metric to use, one of "euclidean" or "correlation"
}
  \item{rm.na}{
    should NA values be imputed?
}
  \item{rm.nan}{
    should NaN values be imputed?
}
  \item{rm.inf}{
    should Inf values be imputed?
}
}
\details{
  Uses the K-nearest neighbor algorithm, as described in Troyanskaya et
  al., 2001, to impute missing values in a data matrix.  Elements are
  imputed row-wise, so that neighbors are
  selected based on the rows which are closest in distance to the row
  with missing values.  There are two
  choices for a distance metric, either Euclidean (the default) or a
  correlation 'metric'.  If the latter is selected, matrix values are
  first row-normalized to mean zero and standard deviation one to select
  neighbors.  Values are 'un'-normalized by applying the inverse
  transformation prior to returning the imputed data matrix. 
}
\value{
  A data matrix with missing values imputed.
}
\references{

  O. Troyanskaya, M. Cantor, G. Sherlock, P. Brown, T. Hastie,
  R. Tibshirani, D. Botstein, and R. B. Altman.
  Missing value estimation methods for dna microarrays.
  Bioinformatics, 17(6):520-5, 2001.

  G.N. Brock, J.R. Shaffer, R.E. Blakesley, M.J. Lotz, and G.C. Tseng.
  Which missing value imputation method to use in expression profiles: a
  comparative study and two selection schemes.
  BMC Bioinformatics, 9:12, 2008.
  
}
\author{
 Guy Brock 
}

\seealso{
  See the package vignette for illustration on usage. 
}
\examples{

## generate some fake data and impute MVs
set.seed(101)
mat <- matrix(rnorm(500), nrow=100, ncol=5)
idx.mv <- sample(1:length(mat), 50, replace=FALSE)
mat[idx.mv] <- NA
imputed <- imputeKNN(mat)

}
\keyword{ manip }

