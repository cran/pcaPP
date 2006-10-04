plotcov <- function (cov1, cov2, method1, labels1, method2, labels2, ndigits = 4, ...)
{

	if (class (cov1) == "matrix")
		cm1 = cov1
	else if (is.null (cov1$cov))
		stop ("No appropriate covariance structure specified")
	else
	{
		cm1 = cov1$cov
		if (!is.null (cov1$method))
			method1 = cov1$method
	}

	if (missing (method1)) #|| is.null(method1))
		method1 = "Method 1"

	if (ncol (cm1) != nrow (cm1))
		stop ("Supplied covariance structure has to be quadratic!")

	if (missing (labels1))
		if (is.null (dimnames (cm1)))
			labels1 = paste ("V", 1:ncol (cm1), sep = "")
		else
			labels1 = dimnames (cm1)[[1]]

	if (missing (cov2))
	{	# only plotting cov1
		plotsingle (cm1, labels1, method1, ndigits, ...)
	}
	else
	{
		if (class (cov2) == "matrix")
			cm2 = cov2
		else if (is.null (cov2$cov))
			stop ("No appropriate covariance structure specified")
		else
		{
			cm2 = cov2$cov
			if (!is.null (cov2$method))
				method2 = cov2$method
		}
		if (missing (method2))# || is.null(method2))
			method2 = "Method 2"

		if (ncol (cm2) != nrow (cm2))
			stop ("Supplied covariance structure has to be quadratic!")

		if (missing (labels2))
			if (is.null (dimnames (cm2)))
				labels2 = paste ("V", 1:ncol (cm2), sep = "")
			else
				labels2 = dimnames (cm2)[[1]]

		plotcomp (cm1, cm2, labels1, labels2, method1, method2, ndigits, col = "#0096FF", ...)
	}
	invisible()
}
