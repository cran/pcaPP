

.l1median_NLM_Hess <-
function (X, maxit = 200, tol = 10^-8, trace = 0, m.init = apply(X, 2, median), msg = 8, method = 1, GFlag = 1, HFlag = 1, Exp = 1, Digits = 6, ...)
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

	ret = .C ("l1median_NLM_Hess", PACKAGE="pcaPP", NAOK = TRUE, 
		npar = as.integer (c(dim (X), maxit, 0, method, 0, msg, trace, GFlag, HFlag, Exp, Digits)),
		dpar = as.double (c(tol, 0)),
		as.double (X),
		med = as.double (m.init)
		)

	if (trace >= 1)
		cat ("l1median returned", ret$npa, ret$dpar, "\r\n") ;

		if (ret$npar[7])
			stop (paste ("nlm optimization returned error code", ret$npar[7]))

	return (list (par = ret$med, value = ret$dpar[2], code = ret$npar[4], iterations = ret$npar [3], time = ret$npar[6]))
}
