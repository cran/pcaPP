covPC = function (x, k = dim(x$loadings)[2], method)
{
	if (!any(class(x) == "princomp"))
		stop ("Invalid parameter \x22k\x22: Data type princomp expected!")
	if (length (x$sdev) != length (x$center))
		warning ("Calculating a rank", length (x$sdev), "- covariance matrix")

	ret = list()
	k = min (ncol (x$loadings), k)

	ret$cov = x$loadings[,1:k] %*% diag (x$sdev[1:k]^2) %*% t(x$loadings[,1:k])
	ret$center = x$center
	if (missing (method))
		ret$method = "Covariance estimation based on PCs"
	else
		ret$method = method
	return (ret)
} 
