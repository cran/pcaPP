.centr <- function(X,m) t (t(X) - m)


.mymean <- function (X)
{
	if (any (!is.na (X)))
		mean (X[!is.na (X)], trim = 0.4)
	else
		1
}

#.scale <- function (X)
#{
#	med <- apply (X, 2, median)
#	dev <- abs (.centr (X, med))
#	dev[dev == 0] <- NA	
#	apply(dev, 2, .mymean)
#}

# .First.lib <-
# function(lib,pkg)
# {	##	2do: delete, using namespaces now; this function is useless.
#    library.dynam("pcaPP",pkg,lib)
#    library(mvtnorm)
#    cat("pcaPP 0.1-1 loaded\n")
# }

.ParseControlStructure = function (control, arguments)
{
	if (!is.list(control))
		stop ("Invalid argument type: control structure must be of type list")

	if (missing(arguments))
		arguments = attributes(control)$names

	for (curname in arguments)
	{
		if (!is.null (control[[curname]]))	##	 if this argument is provided in the control structure:
			eval (parse (text = paste ("eval.parent (substitute(", curname, "<- control$", curname, "), n = 3)", sep = "")))
	}
}

.ParseDevString = function (method)
{
	if (method[1] == "mad")  return (0)
	if (method[1] == "sd")  return (1)
	if (method[1] == "Qn" | method[1] == "qn" )  return (2)
	return (1)
}

.DataPostProc <- function (DataObj, obj, loadings, scores, cl, bScores)
{
	if (bScores)
	{
		idx <- order (obj, decreasing = TRUE)
		scores <- as.matrix (scores [,idx])
		obj <- obj [idx]
		loadings <- as.matrix (loadings [,idx])
	}
	ret <- list()

   ##loadings
	{
		c <- ncol (loadings)
		r <- nrow (loadings)
		ret$loadings <- loadings

		ret$loadings <- .loadSgnU (ret$loadings)

		if (is.null (dimnames (DataObj$x)[[2]]))
			dimnames (ret$loadings) <- list (paste (rep ("V", r), 1:r, sep = ""), paste (rep ("Comp.", c), 1:c, sep = ""))
		else
			dimnames (ret$loadings) <- list (dimnames (DataObj$x)[[2]], paste (rep ("Comp.", c), 1:c, sep = ""))

		class (ret$loadings) <- "loadings"
	}

   ##sdev
	ret$sdev <- as.numeric (obj)
	names (ret$sdev) <- dimnames (ret$loadings)[[2]]

   ##center
	ret$center <- DataObj$center
   ##scale
	ret$scale <- DataObj$scale
   ##n.obs
	ret$n.obs <- nrow (DataObj$x)

   ##scores
	if (bScores)
	{
		ret$scores <- scores
		dimnames (ret$scores) <- list (1:nrow (scores), dimnames (ret$loadings)[[2]]) ;
	}
	else
		ret$scores <- NULL

	ret$call <- cl

	class (ret) <- c ("pcaPP", "princomp")
	return (ret)
}

.Conv2Matrix <- function (x, sx = substitute (x))
{
	if(is.matrix(x))
		return (x)
	if(is.data.frame(x))
		return (data.matrix(x))
	return (matrix(x, nrow = length(x), ncol = 1, dimnames = list(names(x), deparse(sx))))
}


.colMedians <- function (x)
{
	if (is.null (dim (x)))
		return (median (x))
	apply (x, 2, median)
}

.GetFunctionName <- function (f, ...)
{
	form <- formals (f)
	if (!is.null (form$fa.Name))
		return (eval (form$fa.Name))

	if (!is.null (form$NAME))
		return (f (NAME = TRUE, ...))

	if (!is.null (attributes (f)$NAME))
		return (attributes (f)$NAME)

	return (NULL)
}

.flush.cat <- function (...)
{
	cat (...)
	flush.console ()	
}
