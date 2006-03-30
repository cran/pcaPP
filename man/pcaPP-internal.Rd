\name{pcaPP-internal}
\title{Internal pcaPP objects}
\alias{.First.lib}
\alias{DataPostProc}
\alias{EvalScaleFunction}
\alias{ParseControlStructure}
\alias{ParseDevString}
\usage{
DataPostProc(DataObj, obj, loadings, scores, cl, bScores)
EvalScaleFunction(method)
ParseDevString(method)
ParseControlStructure(control, arguments)
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
}
\keyword{internal}
