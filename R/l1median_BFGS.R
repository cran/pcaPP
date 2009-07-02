l1median_BFGS <-
function (X, maxit = 200, tol = 10^-8, trace = 0, m.init = apply(X, 2, median), REPORT = 10, ...)
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
		
		
	ret = .C ("l1median_BFGS", PACKAGE="pcaPP", NAOK = TRUE, 
		par = as.integer (c(dim (X), maxit, trace, REPORT)),
		npar.out = integer (4),
		dpar = as.double (c(-Inf, tol)),
		dpar.out = double (1),
		as.double (X),
		#as.double (pscale),
		med = as.double (m.init)#double (ncol(X))
		)

	if (trace >= 1)
		cat ("l1median returned", ret$npa, ret$dpar, "\r\n") ;

	return (list (par = ret$med, value = ret$dpar.out[1], code = ret$npar.out [1], iterations = ret$npar.out [2], iterations_gr = ret$npar.out [3], time = ret$npar.out[4]))
}

