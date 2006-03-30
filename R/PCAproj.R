"PCAproj" <-
function (x, k = 2, method = c("sd", "mad", "qn"), CalcMethod = c("eachobs", "lincomb", "sphere"), 
          nmax = 1000, update = TRUE, scores = TRUE, maxit = 5, maxhalf = 5, control, ...)
{
	if (!missing (control))
		ParseControlStructure (control, c("k", "method", "CalcMethod", "nmax", "update", "scores", "maxit", "maxhalf"))

	if (!any (CalcMethod == c("eachobs", "lincomb", "sphere")))
		stop (paste ("Unknown calcmethod:", CalcMethod))

	k = k[1]
	method = method[1]
	CalcMethod = CalcMethod[1]
	nmax = nmax[1]
	update = as.logical (update[1])
	scores = as.logical (scores[1])
	maxit = maxit[1]
	maxhalf = maxhalf[1]

	n = nrow (x)
	p = ncol (x)

	if( k > min(n,p))
		stop ('k too large')

        if(p > n)
        {
                svdx = svd(t(x))
                x = svdx$v %*% diag(svdx$d)
                pold=p
                p=n
        } else
                pold=p


	DataObj = ScaleAdvR (x, ...)
	y = DataObj$x

#	m = l1median(x)
#	y = t(t(x) - m)

	if (scores)
		scoresize = n*k
	else
		scoresize = 0

	if (CalcMethod == "lincomb")
	{
		update = FALSE
		if (nmax > n)
		{
			aux = matrix (runif ((nmax-n) * n), nrow = nmax-n)
			y = rbind (y, t(t(aux %*% as.matrix(x)) - DataObj$center))
		}
	}
	else if (CalcMethod == "sphere")
	{
		update = FALSE
		if(nmax>n)
			#y[(n+1):nmax,] =	rmvnorm(nmax-n, rep(0,p), diag (p))
			y = rbind (y,		rmvnorm(nmax-n, rep(0,p), diag (p)))
	}

	nn = nrow (y)

	if (update)
		ret = .C ("rpcnup", as.double (y), as.integer (c(nn, p, k, ParseDevString (method), scores, maxit, maxhalf, n)),
			scores = double (scoresize),
			loadings = double (p * k),
			lambda = double (k),
			PACKAGE = "pcaPP")
	else
		ret = .C ("rpcn",   as.double (y), as.integer (c(nn, p, k, ParseDevString (method), scores, n)),
			scores = double (scoresize),
			loadings = double (p * k),
			lambda = double (k),
			PACKAGE = "pcaPP")


       if(pold>n)
            veig = svdx$u %*% veig

	if (scores)
		DataPostProc (DataObj, ret$lambda, matrix (ret$loadings, ncol = k), matrix (ret$scores, ncol = k) , match.call(), scores)
	else
		DataPostProc (DataObj, ret$lambda, matrix (ret$loadings, ncol = k), NULL, match.call(), scores)
}

