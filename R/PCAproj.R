
PCAproj <- function (x, k = 2, method = c ("mad", "sd", "qn"), CalcMethod = c("eachobs", "lincomb", "sphere"), nmax = 1000, update = TRUE, scores = TRUE, maxit = 5, maxhalf = 5, scale = NULL, center = l1median_NLM, zero.tol = 1e-16, control)
{
	if (!missing (control))
		###ParseControlStructure (control, c("k", "method", "CalcMethod", "nmax", "update", "scores", "maxit", "maxhalf"))
		.ParseControlStructure (control, c("k", "method", "CalcMethod", "nmax", "update", "scores", "maxit", "maxhalf", "center", "scale", "zero.tol"))

	method <- .getScaleMethod (method)

	CalcMethod <- match.arg (CalcMethod[1], c("eachobs", "lincomb", "sphere"))

	x <- .Conv2Matrix (x, substitute (x))

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

	DataObj = ScaleAdv (x, scale = scale, center = center)

    if (pold > n) # center and scale must have original data dimension:
	{
		DataObj$center <- as.vector(svdx$u%*%DataObj$center)
		DataObj$scale <- ScaleAdv(x%*%t(svdx$u),center=NULL,scale=scale)$scale
	}

	y = DataObj$x

#	m = l1median(x)
#	y = t(t(x) - m)

	if (scores)
		scoresize <- n * k
	else
		scoresize <- 0

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
		if(nmax  >n)
			#y[(n+1):nmax,] =	rmvnorm(nmax-n, rep(0,p), diag (p))
			y = rbind (y,		rmvnorm(nmax-n, rep(0,p), diag (p)))
	}

	nn = nrow (y)

	if (update)
		ret.C = .C ("pcaProj_up", PACKAGE="pcaPP"
				, as.integer (c(nn, p, n, k, method, scores, maxit, maxhalf))
				, as.double (zero.tol)
				, as.double (y)
				, scores = double (scoresize)
				, loadings = double (p * k)
				, lambda = double (k)
			)
	else
		ret.C = .C ("pcaProj",  PACKAGE="pcaPP"
				, as.integer (c(nn, p, n, k, method, scores))
				, as.double (zero.tol)
				, as.double (y)
				, scores = double (scoresize)
				, loadings = double (p * k)
				, lambda = double (k)
			)


	veig = matrix (ret.C$loadings, ncol = k)

   if(pold>n)
		veig = svdx$u %*% veig

	if (scores)
		.DataPostProc (DataObj, ret.C$lambda, veig, matrix (ret.C$scores, ncol = k) , match.call(), scores)
	else
		.DataPostProc (DataObj, ret.C$lambda, veig, NULL, match.call(), scores)
}
