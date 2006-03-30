\name{l1median}
\alias{l1median}
\title{ Multivariate L1 Median }
\description{
  Computes the multivariate L1 median (also called spatial median) of a data matrix.  
}
\usage{
l1median(X, MaxStep = 200, ItTol = 10^-8)
}
\arguments{
  \item{X}{ a matrix containing the values whose multivariate L1 median is to be computed.}
  \item{MaxStep}{Maximal number of step-halvings.}
  \item{ItTol}{Tolerance for convergence of the algorithm.}
}
\value{
  returns the vector of the coordinates of the L1 median.
}
\references{ C. Croux, P. Filzmoser, M. Oliveira (2004)
  \emph{Projection-pursuit Estimators for Robust Principal Component Analysis,
        Technical Report TS-04-4, Vienna University of Technology, Austria}
}
\author{Heinrich Fritz, Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>}

\seealso{ \code{\link{median}}}
\examples{
  l1median(rnorm(100)) # this returns the median of the sample

  # multivariate data with outliers
  library(mvtnorm)
  x <- rbind(rmvnorm(200, rep(0, 4), diag(c(1, 1, 2, 2))), 
             rmvnorm( 50, rep(3, 4), diag(rep(2, 4))))
  l1median(x)
  # compare with coordinate-wise median:
  apply(x,2,median)
}
\keyword{ multivariate }
\keyword{ robust }
