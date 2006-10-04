covPCAgrid = function (x, control)
{
	pcs = PCAgrid (x, k = ncol(x), control = control)

	ret = list()
	ret$cov = pcs$loadings %*% diag (pcs$sdev^2) %*% t(pcs$loadings)
	ret$center = pcs$center
	ret$method = "Robust cov - estimation based on PCs (grid mode)"

	if (!missing (control) && !is.null (control$method))
		ret$method = paste ("Robust cov - estimation based on PCs (grid mode -", control$method, ")", sep = "")
	else
		ret$method = "Robust cov - estimation based on PCs (grid mode)"

	return (ret)
}
