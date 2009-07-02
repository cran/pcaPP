l1median_HoCr <-
function (X, maxit = 200, tol = 10^-8, trace = 0, m.init = apply(X, 2, median), ...)
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

	ret = .C ("l1median_HoCr", PACKAGE="pcaPP", 
		npar = as.integer (c(dim (X), maxit, trace)),
		npar.out = integer (4),
		dpar = as.double (c (tol)),
		dpar.out = double (1),
		as.double (X),
		med = as.double (m.init)
		)

	if (ret$npar[3] == 2)
		warning ("step halving failed")# after ... steps")

	return (list (par = ret$med, value = ret$dpar.out[1],
		 code = ret$npar.out[3], iterations = ret$npar.out [1] + 1, stephalving = ret$npar.out [2], time = ret$npar.out[4]))
}
