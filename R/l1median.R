l1median = function (X, MaxStep = 200, ItTol = 10^-8, trace = 0)
{
	if (trace >= 0)
		warning ("This function (pcaPP::l1median)is outdated.\r\nFor better performance try any of pcaPP::l1median_* instead. Preferably pcaPP::l1median_HoCrJo.\r\nOtherwise use (trace = -1) for suppressing this warning. ")

	if (class (X) != "matrix")
	{
		if (class (X) == "data.frame")
			X = as.matrix(X)
		else
			X = matrix(X, ncol = 1)
	}

	ret = .C ("l1median", PACKAGE="pcaPP",
		as.double (X),
		as.integer (nrow(X)),
		as.integer (ncol(X)),
		med = double (ncol(X)),
		ret = integer(1),
		as.integer (MaxStep),
		as.double (ItTol)
		)

	if (ret$ret != 0)
		return (ret$med)
	stop("iteration failed")
}
