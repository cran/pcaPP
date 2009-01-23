\name{qn}
\alias{qn}
\title{ scale estimation using the robust Qn estimator }
\description{
  Returns a scale estimation as calculated by the (robust) Qn estimator.
}
\usage{
qn(x)
}
\arguments{
  \item{x}{ a vector of data }
}
\details{
  The Qn estimator computes the first quartile of the pairwise absolute 
  differences of all data values.
}
\value{
  The estimated scale of the data.
}
\references{ 
  P.J. Rousseeuw, C. Croux (1993)
  Alternatives to the Median Absolute Deviation,
  \emph{JASA}, \strong{88}, 1273-1283.
}
\author{ Heinrich Fritz, Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>}

\seealso{ \code{\link{mad}} }
\examples{
  # data with outliers
  x <- c(rnorm(100), rnorm(10, 10))
  qn(x)
}
\keyword{multivariate}
\keyword{robust}
