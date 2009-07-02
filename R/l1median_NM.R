l1median_NM <-
function (X, maxit = 200, tol = 10^-8, trace = 0, m.init = apply(X, 2, median), pscale = apply(abs(.centr(X, m.init)), 2, mean, trim = 0.40))
{
#	centr <- function(X,m) X - rep(m, each = nrow (X))

	if (class (X) != "matrix")
	{
		if (class (X) == "data.frame")
			X = as.matrix(X)
		else
			X = matrix(X, ncol = 1)
	}
	
	ret = .C ("l1median_NM", PACKAGE="pcaPP", NAOK = TRUE, 
		npar = as.integer (c(dim (X), maxit, 0, 0, 0)),
		dpar = as.double (c(-Inf, tol, 0, 1, 0.5, 2)),
		as.double (X),
		as.double (pscale),
		med = as.double (m.init)#double (ncol(X))
		)

	if (trace >= 1)
		cat ("l1median returned", ret$npa, ret$dpar, "\r\n") ;

	return (list (par = ret$med, value = ret$dpar[3], code = ret$npar[4], iterations = ret$npar [6]))
}

