.centr <-
function(X,m) t (t(X) - m)


.mymean <-
function (X)
{
	if (any (!is.na (X)))
		mean (X[!is.na (X)], trim = 0.4)
	else
		1
}

.scale <-
function (X)
{
	med <- apply (X, 2, median)
	dev <- abs (.centr (X, med))
	dev[dev == 0] <- NA	
	apply(dev, 2, .mymean)
	
	#pscale <- apply(dev, 2, mean, trim = 0.40, rm.na = T)


#	pscale [is.nan (pscale )] = 1		##	if any variables variance is zero it's scale value is set to 1 (NaN resulting from mean (c (NA, Na, Na, ... , NA), nr.rm = T))
#	pscale
}

 .First.lib <-
 function(lib,pkg)
 {
    library.dynam("pcaPP",pkg,lib)
    library(mvtnorm)
    cat("pcaPP 0.1-1 loaded\n")
 }

