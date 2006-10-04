\name{pcaPP-internal}
\title{Internal pcaPP objects}
\alias{.First.lib}
\alias{DataPostProc}
\alias{EvalScaleFunction}
\alias{ParseControlStructure}
\alias{ParseDevString}
\alias{doEllipses}
\alias{plotcomp}
\alias{plotsingle}
\usage{
DataPostProc(DataObj, obj, loadings, scores, cl, bScores)
EvalScaleFunction(method)
ParseDevString(method)
ParseControlStructure(control, arguments)
doEllipses(acov, pos, ...)
plotcomp(cor1, cor2, labels1, labels2, method1, method2, ndigits = 4, ...)
plotsingle(cm, labels, method, ndigits, ...)
 
}
\description{Internal pcaPP objects.}
\details{These are not to be called by the user.

DataPostProc creates a pcaPP object using the given parameters. It is called
to format the output to the correct shape.

EvalScaleFunction and ParseDevString are used to handle parameters that 
should be functions for calculating the (robust) standard deviation.

ParseControlStructure handles the control object parameter used in 
\code{\link{covPCAproj}}, \code{\link{covPCAgrid}}, \code{\link{PCAproj}}
and \code{\link{PCAgrid}}.

doEllipses, plotcomp, and plotsingle are used in \code{\link{plotcov}}
to visually compare two covariance matrix estimations.
}
\keyword{internal}
