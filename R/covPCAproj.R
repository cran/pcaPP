"covPCAproj" <-
function (x, control, ...)
{	## es ist nicht möglich die covPC - funktion einfach auf das princomp - objekt anzuwenden, da covrob, damit umgehen können  muss, und dafür muss die datenmatrix übergeben werden!
	pcs = PCAproj (x, k = ncol(x), control = control, ...)

	ret = list()
	ret$cov = pcs$loadings %*% diag (pcs$sdev^2) %*% t(pcs$loadings)
	ret$center = pcs$center
	ret$method = "Robust cov - estimation based on PCs (projection mode)"
	return (ret)
}

