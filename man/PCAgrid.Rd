\name{PCAgrid}
\alias{PCAgrid}
\title{ Robust Principal Components using the Grid search algorithm }
\description{
  Computes a desired number of (robust) principal components using the grid search 
  algorithm in the plane. The global optimum of the objective function is 
  searched in planes, not in the p-dimensional space, using regular grids in 
  these planes.
}
\usage{
PCAgrid(x, k = 2, method = c("mad", "sd", "qn"), maxiter = 10, splitcircle = 10, 
scores = TRUE, anglehalving = TRUE, fact2dim = 10, scale = NULL, center = l1median, control)
}
\arguments{
  \item{x}{ a numeric matrix or data frame which provides the data for the 
  principal components analysis. }
  \item{k}{ desired number of components to compute }
  \item{method}{ scale estimator used to detect the direction with the largest
  variance. Possible values are \code{"sd"}, \code{"mad"} and \code{"qn"}, the 
  latter can be called \code{"Qn"} too. \code{"mad"} is the default value.}
  \item{maxiter}{ maximum number of iterations. }
  \item{splitcircle}{ the number of directions in which the algorithm should
  search for the largest variance. The direction with the largest variance
  is searched for in the directions defined by a number of equally spaced points
  on the unit circle. This argument determines, how many such points are used to
  split the unit circle. }
  \item{scores}{ a logical value indicating whether the scores of the
  principal component should be calculated. }
  \item{anglehalving}{ boolean stating whether angle halving is to be used or not. 
  Angle halving will usually improve the solution quite a lot.}
  \item{fact2dim}{ an integer that is multiplied to splitcircle if x is only 
  two-dimensional. In higher dimensions, fewer search directions are needed to allow
  for faster computation. In two dimensions, more search directions are required to 
  grant higher precision. \code{fact2dim} is used to take account of this.}
  \item{scale}{ this argument indicates how the data is to be rescaled. It
  can be a function like \code{\link{sd}} or \code{\link{mad}} or a vector
  of length \code{ncol(x)} containing the scale value of each column. }
  \item{center}{ this argument indicates how the data is to be centered. It
  can be a function like \code{\link{mean}} or \code{\link{median}} or a vector
  of length \code{ncol(x)} containing the center value of each column. }
  \item{control}{ a list whose elements must be the same as (or a subset of) 
  the parameters above. If the control object is supplied, the parameters from 
  it will be used and any other given parameters are overridden. }
}
\details{
  Angle halving is an extension of the original algorithm. In the original 
  algorithm, the search directions are determined by a number of points on the
  unit circle in the interval [-pi/2 ; pi/2). Angle halving means this angle is
  halved in each iteration, eg. for the first approximation, the above mentioned
  angle is used, for the second approximation, the angle is halved to 
  [-pi/4 ; pi/4) and so on. This usually gives better results with less 
  iterations needed.

  Similar to the function \code{\link{princomp}}, there is a \code{print} method 
  for the these objects that prints the results in a nice format and the 
  \code{plot} method produces a scree plot (\code{\link{screeplot}}). There is 
  also a \code{biplot} method.
}
\value{
  The function returns an object of class \code{"princomp"}, i.e. a list 
  similar to the output of the function \code{\link{princomp}}.
  \item{sdev}{the (robust) standard deviations of the principal components.}
  \item{loadings}{the matrix of variable loadings (i.e., a matrix whose columns
  contain the eigenvectors).  This is of class \code{"loadings"}:
  see \code{\link{loadings}} for its \code{\link{print}} method.}
  \item{center}{the means that were subtracted.}
  \item{scale}{the scalings applied to each variable.}
  \item{n.obs}{the number of observations.}
  \item{scores}{if \code{scores = TRUE}, the scores of the supplied data on the
  principal components.}
  \item{call}{the matched call.}
}
\references{
  C. Croux, P. Filzmoser, M. Oliveira, (2007).
  Algorithms for Projection-Pursuit Robust Principal Component Analysis,
  \emph{Chemometrics and Intelligent Laboratory Systems}. To appear.
}
\author{Heinrich Fritz, Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>}
\seealso{ \code{\link{PCAproj}}, \code{\link{ScaleAdv}}, \code{\link{princomp}} }
\examples{
  # multivariate data with outliers
  library(mvtnorm)
  x <- rbind(rmvnorm(200, rep(0, 6), diag(c(5, rep(1,5)))),
             rmvnorm( 15, c(0, rep(20, 5)), diag(rep(1, 6))))
  # Here we calculate the principal components with PCAgrid
  pc <- PCAgrid(x)
  # we could draw a biplot too:
  biplot(pc)
  # now we want to compare the results with the non-robust principal components
  pc <- princomp(x)
  # again, a biplot for comparison:
  biplot(pc)
}
\keyword{multivariate}
\keyword{robust}
