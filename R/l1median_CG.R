l1median_CG <-
function (X, maxit = 200, tol = 10^-8, trace = 0, m.init = apply(X, 2, median), pscale = apply(abs(.centr(X, m.init)), 2, mean, trim = 0.40))#, type = 1)
{
	## 2do: check influence of type on the result / runtime
	##			all types SEEM to give the same result, type 1 SEEMS to be the fastest!
	type = 1	##	 type shall be not an option anymore

#	centr <- function(X,m) X - rep(m, each = nrow (X))

	if (class (X) != "matrix")
	{
		if (class (X) == "data.frame")
			X = as.matrix(X)
		else
			X = matrix(X, ncol = 1)
	}

	if (type < 1 || type > 3)
		stop ("parameter type MUST be either 1, 2 or 3")

		ret = .C ("l1median_CG", PACKAGE="pcaPP", NAOK = TRUE, 
		par = as.integer (c(dim (X), maxit, trace, type)),
		npar.out = integer (3),
		dpar = as.double (c(-Inf, tol)),
		dpar.out = double (1),
		as.double (X),
		as.double (pscale),
		med = as.double (m.init)#double (ncol(X))
		)

	if (trace >= 1)
		cat ("l1median returned", ret$npa, ret$dpar, "\r\n") ;

	return (list (par = ret$med, value = ret$dpar.out[1], code = ret$npar.out [1], iterations = ret$npar.out [2], iterations_gr = ret$npar.out [3]))
}

