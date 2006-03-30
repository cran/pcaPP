#include "math.h"
#include "R.h"

#include "fastpca.h"
#include "rsubst.h"


/*

  "scale.tau"<-function(y, center = median(y), weights = NULL, init.scale = NULL, tuning = 1.95, na.rm = F)
{
	if(!missing(weights))
	{
		if(length(weights) != length(y))
			stop("length of weights must be the same as length of y"
				)
		if(any(is.na(weights)))
			stop("missing values not allowed in weights")
		if(any(weights < 0))
			stop("negative weights not allowed")
	}
	else
		weights <- rep(1, length(y))
	nas <- is.na(y)
	if(any(nas))
		if(na.rm)
		{
			y <- y[!nas]
			weights <- weights[!nas]
		}
		else return(NA)
	y <- y - center
	if(is.null(init.scale))
		init.scale <- mad(y, center = 0, low = T)
	if(init.scale <= 0)
		return(0)
	y <- y/init.scale
	if(tuning == 1.95) return(1.048 * init.scale * sqrt(sum(weights * pmin(y^2, 3.8025))/sum(weights)))	#
		# compute constant for arbitrary tuning parameter
	ff <- function(x)
	pmin(x^2, tuning^2) * dnorm(x)
	xs <- seq(-10, 10, length = 1000)
	h <- xs[2] - xs[1]
	assign("tuning", tuning, frame = 1)
	consistency <- h * sum(c(17/48, 59/48, 43/48, 49/48, rep(1, length = 992), 49/48, 43/48, 59/48, 17/48) * ff(xs))
	init.scale/sqrt(consistency) * sqrt(mean(pmin(y^2, tuning^2)))
}

*/

	extern bool scale_tau (double *pdData, int *pnSize, double dCenter, double *pdWeights, double dInit_Scale, double dTuning, int bNa_rm)
	{
/*		if (pdWeights)
		{
			if (VectorMin (pdWeights) < 0)
				return FALSE ;
		}
			//	na gibts nicht..
		double *pdDataC = NULL ;
		
		if (dCenter)
		{
			pdDataC = new double [*pnSize] ;
			matcpy (pdDataC, pdData, *pnSize) ;
			VectorAdd (pdDataC, pnSize, -dCenter) ;
			pdData = pdDataC ;
		}

		

		if (pdDataC)
			delete [] pdDataC ;*/
		return TRUE ;
	}

/*

  "scale.a"<-
function(y, center = median(y), weights = NULL, init.scale = NULL, tuning = 3.85, na.rm = F)
{
	if(!missing(weights)) {
		if(length(weights) != length(y))
			stop("length of weights must be the same as length of y"
				)
		if(any(is.na(weights)))
			stop("missing values not allowed in weights")
		if(any(weights < 0))
			stop("negative weights not allowed")
	}
	else weights <- rep(1, length(y))
	nas <- is.na(y)
	if(any(nas))
		if(na.rm) {
			y <- y[!nas]
			weights <- weights[!nas]
		}
		else return(NA)
	y <- y - center
	if(is.null(init.scale))
		init.scale <- mad(y, center = 0, low = T)
	if(init.scale <= 0)
		return(0)
	y <- y/init.scale
	if(tuning == 3.85) {
		num <- ifelse(abs(y) > 3.85, 0, y * y * (14.8225 - y * y)^4)
		den <- ifelse(abs(y) > 3.85, 0, 219.7065 - 88.935 * y * y + 5 * 
			y^4)
		return((0.9471 * sqrt(sum(weights)) * init.scale * sqrt(sum(
			weights * num)))/abs(sum(weights * den)))
	}
#
# compute constants for arbitrary tuning parameter
	fh <- function(x)
	x * x * (tuning * tuning - x * x)^4 * dnorm(x)
	xs <- seq( - tuning, tuning, length = 600)
	h <- xs[2] - xs[1]
	assign("tuning", tuning, frame = 1)
	h.int2 <- h * sum(c(17/48, 59/48, 43/48, 49/48, rep(1, length = 592), 
		49/48, 43/48, 59/48, 17/48) * fh(xs))
	fg <- function(x)
	(tuning^4 - 6 * tuning * tuning * x * x + 5 * x^4) * dnorm(x)
	g.int <- h * sum(c(17/48, 59/48, 43/48, 49/48, rep(1, length = 592), 49/
		48, 43/48, 59/48, 17/48) * fg(xs))
	num <- ifelse(abs(y) > tuning, 0, y * y * (tuning * tuning - y * y)^4)
	den <- ifelse(abs(y) > tuning, 0, tuning^4 - 6 * tuning * tuning * y * 
		y + 5 * y^4)
	(g.int/sqrt(h.int2) * sqrt(sum(weights)) * init.scale * sqrt(sum(
		weights * num)))/abs(sum(weights * den))
}

*/  
