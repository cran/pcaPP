l1median_NM <-
#function (X, maxit = 200, tol = 10^-8, trace = 0, m.init = apply(X, 2, median), alpha = 1, beta = 0.5, gamma = 2, ...)
function (X, maxit = 200, tol = 10^-8, trace = 0, m.init = apply(X, 2, median), ...)
{
	alpha = 1
	beta = 0.5
	gamma = 2
	
	if (class (X) != "matrix")
	{
		if (class (X) == "data.frame")
			X = as.matrix(X)
		else
			X = matrix(X, ncol = 1)
	}
	
	if (length (m.init) != ncol (X))
		stop (paste ("length of vector m.init (=", length (m.init), ") does not match the number of columns of data object X (=", ncol (X),")", sep = ""))
		
	ret = .C ("l1median_NM", PACKAGE="pcaPP", NAOK = TRUE, 
		npar = as.integer (c(dim (X), maxit, 0, 0, 0, 0, trace)),
		dpar = as.double (c(-Inf, tol, 0, alpha, beta, gamma)),
		as.double (X),
		#as.double (pscale),
		med = as.double (m.init)#double (ncol(X))
		)

	if (trace >= 1)
		cat ("l1median returned", ret$npa, ret$dpar, "\r\n") ;

	return (list (par = ret$med, value = ret$dpar[3], code = ret$npar[4], iterations = ret$npar [6], time = ret$npar[7]))
}

