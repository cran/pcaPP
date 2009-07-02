l1median_VaZh <- 
function (X, maxit = 200, tol = 10^-8, zero.tol = 1e-15, trace = 0, m.init = apply (X, 2, median), ...)
{
	if (class (X) != "matrix")
	{
		if (class (X) == "data.frame")
			X = as.matrix(X)
		else
			X = matrix(X, ncol = 1)
	}
	
	if (length (m.init) != ncol (X))
		stop (paste ("length of vector m.init (=", length (m.init), ") does not match the number of columns of data object X (=", ncol (X),")", sep = ""))

	ret = .C ("l1Median_VZ", PACKAGE="pcaPP", 
		npar = as.integer (c(dim (X), maxit, 0, trace, 0)),
		dpar = as.double (c (tol, zero.tol)),
		as.double (X),
		med = as.double (m.init)
		)

	return (list (par = ret$med, value = sum (sqrt (colSums ((t(X) - ret$med)^2))),
		 code = ret$npar[4], iterations = ret$npar [3], time = ret$npar[6]))
}


