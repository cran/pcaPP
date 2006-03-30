"PCAgrid" <-
function (x, k = 2, method = c("sd", "mad", "qn"), maxiter = 10, splitcircle = 10, 
          scores = TRUE, anglehalving = TRUE, fact2dim = 10, control, ...)
{
# grid search for PCA
#
# x ... centered (!) data matrix
# k ... number of components to determine
# m ... number of directions to search ## umbenannt auf "splitcircle"
# maxiter ... maximum number of iterations
# method ... how to estimate the standard deviation
#        should be: "sd", "mad", qn
#

	if (!missing (control))
		ParseControlStructure (control, c("k", "method", "maxiter", "splitcircle", "scores", "anglehalving", "fact2dim"))

	k = k[1]
	method = method[1]
	maxiter = maxiter[1]
	splitcircle = splitcircle[1]
	scores = as.logical (scores[1])
	anglehalving = as.logical (anglehalving[1])
	fact2dim = fact2dim[1]

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


	DataObj = ScaleAdvR (x, ...)

	if (scores)
		res = .C("rpcgrid",
			as.double (DataObj$x),
			as.integer (c(nrow (DataObj$x), ncol(DataObj$x), k, splitcircle, maxiter, ParseDevString (method), anglehalving, scores, fact2dim)),
			obj = double(k),
			loadings = double (ncol(DataObj$x) * k),
			scores = double (nrow(DataObj$x) * k),
			PACKAGE = "pcaPP"
			)
	else
		res = .C("rpcgrid",
			as.double (DataObj$x),
			as.integer (c(nrow (DataObj$x), ncol(DataObj$x), k, splitcircle, maxiter, ParseDevString (method), anglehalving, scores, fact2dim)),
			obj = double(k),
			loadings = double (ncol(DataObj$x) * k),
			scores = double (1),
			PACKAGE = "pcaPP"
			)
	veig = matrix (res$loadings, ncol = k)

	if(pold>n)
		veig = svdx$u %*% veig

	return (DataPostProc (DataObj, res$obj, loadings = veig, scores = matrix (res$scores, ncol = k), match.call(), scores))
}

