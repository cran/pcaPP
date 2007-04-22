PCAgrid = function (x, k = 2, method = c("mad", "sd", "qn"), maxiter = 10, splitcircle = 10, scores = TRUE, anglehalving = TRUE, fact2dim = 10, scale = NULL, center = l1median, control)
{
# grid search for PCA
#
# x ... centered (!) data matrix
# k ... number of components to determine
# m ... number of directions to search ## umbenannt auf "splitcircle"
# maxiter ... maximum number of iterations
# method ... how to estimate the standard deviation
#        should be: "sd", "mad"
#

	if (!missing (control))
		ParseControlStructure (control, c("k", "method", "maxiter", "splitcircle", "scores", "anglehalving", "fact2dim"))
	method = method[1]

	n = nrow (x)
	p = ncol (x)

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


	if (scores)
		res = .C("rpcgrid", PACKAGE="pcaPP",
			as.double (DataObj$x),
			as.integer (c(nrow (DataObj$x), ncol(DataObj$x), k, splitcircle, maxiter, ParseDevString (method), anglehalving, scores, fact2dim)),
			obj = double(k),
			loadings = double (ncol(DataObj$x) * k),
			scores = double (nrow(DataObj$x) * k)
			)
	else
		res = .C("rpcgrid", PACKAGE="pcaPP",
			as.double (DataObj$x),
			as.integer (c(nrow (DataObj$x), ncol(DataObj$x), k, splitcircle, maxiter, ParseDevString (method), anglehalving, scores, fact2dim)),
			obj = double(k),
			loadings = double (ncol(DataObj$x) * k),
			scores = double (1)
			)
	veig = matrix (res$loadings, ncol = k)

	if(pold>n)
		veig = svdx$u %*% veig

	return (DataPostProc (DataObj, res$obj, loadings = veig, scores = matrix (res$scores, ncol = k), match.call(), scores))
}
