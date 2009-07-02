#include "math.h"
#include "R.h"
#include "rsubst.h"
#include "fastpca.h"


#include "R_ext/Applic.h"


//#include "rsubst.c"
/*
	# Function for grid search in a plane
	"gridplane" <- function(y,m=m,method=method)
	{
		# y ... matrix with two columns of data
		# m ... number of direstions to search
		# method ... how to estimate the standard deviation
		#        should be: "std", "mad"

		nangle <- seq(from=0,by=pi/m,length=m)
		alpha <- cbind(cos(nangle),sin(nangle))
		if (method=="std")
			obj <- apply(y %*% t(alpha), 2, sd)
		else if (method=="mad")
			obj <- apply(y %*% t(alpha), 2, mad)
		else
			stop("Method not found")
		ind <- order(obj)
		objmax <- obj[ind[m]]
		alphamax <- alpha[ind[m],]
		list(objmax=objmax,alphamax=alphamax)
	}*/


	extern void gridplane (double *pdY, int *npYSize, int *pm, dev_Type *pfDev, double *pdObjMax, double *pdAlphaMax)
	{
		//	pdY ..... matrix with two columns of data			//	nYSize x 2
		//	nYSize .. RowCount of pdY
		//	m ....... number of direstions to search
		//	bMethod . how to estimate the standard deviation should be: TRUE = "std", FALSE = "mad"

		int nYSize = *npYSize ;
		int m = *pm ;
//		int bMethod = *pbMethod ;

		double *pdAngle = new double [m] ;			//	m x 1
		double *pdAlpha = new double [m * 2] ;		//	m x 2


		int i ;
		for (i = 0; i < m; i++)
		{
			pdAngle [i] = M_PI  * i / m ;
			pdAlpha [i] = cos (pdAngle [i]) ;
			pdAlpha [i + m] = sin (pdAngle [i]) ;
		}

		double *pdtAlpha = new double [m * 2] ;		//	2 x m

		int nTwoCols = 2 ;

		t (pdAlpha,&m, &nTwoCols, pdtAlpha) ;

		double *pdtFatMat = new double [nYSize * m] ;//(double *) malloc (nYSize * m * sizeof(double)) ;
		matmult (pdY, npYSize, &nTwoCols, pdtAlpha, &m, pdtFatMat) ;
		
		double *pdObj = new double [m] ;//(double *) malloc (m * sizeof(double)) ;
/*		if (bMethod)
			colsd (pdtFatMat, npYSize, &m, pdObj) ;
		else
			colmad (pdtFatMat, npYSize, &m, pdObj) ;
*/
		coldev (pdtFatMat, npYSize, &m, pdObj, pfDev) ;

		int nMaxIdx ;
		wheremax (pdObj, &m, &nMaxIdx) ;

		if (pdObjMax)
			*pdObjMax = pdObj [nMaxIdx] ;
		if(pdAlphaMax)
		{
			pdAlphaMax[0] = pdAlpha [nMaxIdx] ;
			pdAlphaMax[1] = pdAlpha [nMaxIdx + m] ;
		}		

		delete [] pdAngle  ;
		delete [] pdAlpha ;

		delete [] pdtAlpha ;
		delete [] pdtFatMat ;

		delete [] pdObj ;
	}
/*
	"gridplane1" <- function(y,a,m=m)
	{
		# y ... matrix with two columns of data
		# m ... number of direstions to search
		#

		astart <- runif(1,0,pi/m)
		nangle <- seq(from=astart,by=pi/m,length=m)
		alpha <- cbind(cos(nangle),sin(nangle))
		alpha1 <- alpha/+(1+2*alpha[,1]*alpha[,2]*a)

		obj <- apply(y %*% t(alpha1), 2, sd)
		ind <- order(obj)
		objmax <- obj[ind[m]]
		alphamax <- alpha1[ind[m],]
		list(objmax=objmax,alphamax=alphamax)		
	}
*/
/*
	void gridplane1 (double *pdY, int *npYSize, int *pm, int *pa, double *pdObjMax, double *pdAlphaMax)
	{
		int m = *pm ;
		double dStart = rand() * M_PI / m ;

		double *pdAngle = (double *) malloc (m * sizeof(double)) ;			//	m x 1

		double *pdAlpha = (double *) malloc (m * 2 * sizeof(double)) ;		//	m x 2

		int i ;
		for (i = 0; i < m; i++)
		{
			pdAngle [i] = M_PI  * i / m ;
			pdAlpha [i] = cos (pdAngle [i]) ;
			pdAlpha [i + m] = sin (pdAngle [i]) ;
		}

		// not implemented!!!!


	}*/
/*
	# Function for grid search in a plane
	"gridplane2" <- function(y,a,m=m,div=1,method=method)
	{
		# y ... matrix with two columns of data
		# m ... number of direstions to search
		# method ... how to estimate the standard deviation
		#        should be: "std", "mad"

		nangle <- seq(from=(-(m-1)/2),							##	Sequence
		to=((m-1)/2))*pi/div/(m-0.9)   # matlab vers
		alpha <- cbind(cos(nangle),sin(nangle))						##	sin, cos
		alpha1 <- alpha/sqrt(1+2*alpha[,1]*alpha[,2]*a)

		if (method=="std")
			obj <- apply(y %*% t(alpha1), 2, sd)					##	transpose, mat. mult., SpaltenSd
		else if (method=="mad")
			obj <- apply(y %*% t(alpha1), 2, mad)					##	transpose, mat. mult., SpaltenMad
		else
			stop("Method not found")

		ind <- order(obj)
		objmax <- obj[ind[m]]			##	number
		alphamax <- alpha1[ind[m],]		##	vector
		list(objmax=objmax,alphamax=alphamax)
	}
*/

	extern void gridplane2 (double *pdY, int *npYSize, double *pa, int *pm, double *pdDiv, dev_Type *pfDev, double *pdObjMax, double *pdAlphaMax)
	{				//	verified
		int m = *pm ;
		double *pdAngle = (double *) malloc (m * sizeof(double)) ;			//	m x 1
		double *pdAlpha = (double *) malloc (m * 2 * sizeof(double)) ;		//	m x 2
		double *pdAlpha1 = (double *) malloc (m * 2 * sizeof(double)) ;			//	m x 1

		int i ;
		double a = *pa ;
		double dCurDiv ;
		for (i = 0; i < m; i++)
		{
			pdAngle [i] = (-(m-1)/2.0 + i) * M_PI / *pdDiv / (m-0.9) ;
			pdAlpha [i] = cos (pdAngle [i]) ;
			pdAlpha [i + m] = sin (pdAngle [i]) ;

			dCurDiv = ::sqrt ((double) (1+2 * pdAlpha [i] * pdAlpha [i + m] * a)) ;

			pdAlpha1 [i] = pdAlpha [i] / dCurDiv ;
			pdAlpha1 [i + m] = pdAlpha [i + m] / dCurDiv ;
		}

		double *pdtAlpha1 = (double *) malloc (m * 2 * sizeof(double)) ;		//	2 x m

		int nTwoCols = 2 ;

		t (pdAlpha1,&m, &nTwoCols, pdtAlpha1) ;

//		printmat (pdAlpha, m, 2) ;
//		printmat (pdrAlpha, 2, m) ;

		double *pdtFatMat = (double *) malloc (*npYSize * m * sizeof(double)) ;
		matmult (pdY, npYSize, &nTwoCols, pdtAlpha1, &m, pdtFatMat) ;

		double *pdObj = (double *) malloc (m * sizeof(double)) ;
//		if (*pnMethod)
//			colsd (pdtFatMat, npYSize, &m, pdObj) ;
//		else
//			colmad (pdtFatMat, npYSize, &m, pdObj) ;

		coldev (pdtFatMat, npYSize, &m, pdObj, pfDev) ;

		free (pdtFatMat) ;

		int nMaxIdx ;
		wheremax (pdObj, &m, &nMaxIdx) ;

		if (pdObjMax)
			*pdObjMax = pdObj [nMaxIdx] ;
		if(pdAlphaMax)
		{
			pdAlphaMax[0] = pdAlpha1 [nMaxIdx] ;
			pdAlphaMax[1] = pdAlpha1 [nMaxIdx + m] ;
		}		

		free (pdObj) ;
		free (pdAngle) ;
		free (pdAlpha) ;
		free (pdAlpha1) ;
		free (pdtAlpha1) ;
	}

	void calcvectorsqrt (double *pdData, int n, double *pdSqrt)
	{
		int i ;
		for (i = n; i > 0; i--)
			*pdSqrt++ = ::sqrt ((double) (*pdData++)) ;
	}

	double calcsumsqrVectorDiff (double *pD1, double *pD2, int nCount)
	{
		double dRet = 0, dCur ; 
		int i ;

		for (i = nCount - 1; i >= 0 ; i--)
		{
			dCur = pD1[i] - pD2[i] ;
			dRet += dCur * dCur ;
		}
		return dRet ;
	}

	double calcsqrtsumsqr (double *pData, int nCount)
	{
		int i ;
		double dSum = 0 ;
		for (i = nCount - 1; i >= 0; i--)
			dSum += pData [i] * pData[i] ;
		return ::sqrt (dSum) ;
	}

	int VectorLargerThan (double *pdData, int n, int *pnTF, double pMargin)
	{
		int nRet = 0 ;
		int i ;
		for (i = n; i > 0; i--)
		{
			//if (*pnTF++ = (*pdData++ > pMargin))
			*pnTF = (*pdData++ > pMargin) ;

			if (*pnTF++)
				nRet ++ ;
		}
		return nRet ;
	}

	double calcrowsumsq (double *pdData, int n, int p, double *pdRowSumSq)
	{
		int i, j ;
		double dCurSum = 0 ;
		double *pCurData ;
		for (i = 0; i < n; i++)
		{
			double 	&dCurSum = pdRowSumSq [i] ;
			dCurSum = 0 ;
			pCurData = pdData + i ;
			for (j = 0; j < p; j++)
			{
				dCurSum += *pCurData * *pCurData ;
				pCurData += n ;
			}
		}
		return dCurSum ;
	}

/*
	extern void pcagridc (double *px, int *pnParams, double *pObj, double *pLoadings, double *pScores)  //#da als character def " "
	{
		int n  = pnParams [0] ;		//	nRows
		int p  = pnParams [1] ;		//	nCols
		
		int k	= pnParams [2] ;
		int m	= pnParams [3] ;
		int maxiter = pnParams [4] ;
		int	method = pnParams [5] ;

		// grid search for PCA
		//
		// x ... centered (!) data matrix						n x p
		// k ... number of components to determine
		// m ... number of directions to search
		// maxiter ... maximum number of iterations
		// method ... how to estimate the standard deviation
		//        should be: "std", "mad"
		//


		int i ;
		for (i = k - 1; i >= 0; i--)
			pObj[i] = R_NaReal ;

		double *z = pScores ;//= (double *) malloc (n * k * sizeof (double)) ;		//	n x k	//	matrix with scores
		fill (pScores, n * k, R_NaReal) ;

		double *l = pLoadings ;//= (double *) malloc (p * k * sizeof (double)) ;		//	p x k	//	matrix with loadings
		fill (pLoadings, p * k, R_NaReal) ;

		if (p == 2)
		{
			int nCurM = m * 10 ;
			gridplane (px, &n, &nCurM, &method, pObj, l) ;

			if (k > 1)
			{
				l[0 + 1 * p] = -l[1] ;
				l[1 + 1 * p] = l[0] ;

				double *pFatMat = NEWDARRAY (n) ;		//	n x 1
				
				int nOneCol = 1 ;
				matmult (px, &n, &p, l + p, &nOneCol, pFatMat) ;

				if (method)
					colsd (pFatMat, &n, &nOneCol, pObj + 1) ;
				else
					colmad (pFatMat, &n, &nOneCol, pObj + 1) ;

				
				free (pFatMat) ;
			}
		}
		else if (p > 2)
		{
			double *pr = NEWDARRAY (p * p);				//	 p x p
			fillDiag (pr, &p, &p) ;

			double *y = NEWDARRAY (n * p) ;				//	n x p
			double *afin = NEWDARRAY (p) ;				//	p x 1


			int nb ;

			double *pdYColSd = NEWDARRAY (p) ;
			int *pnOrder = NEWNARRAY (p) ;

			double *pdAfinBest = NEWDARRAY (p) ;

			double *pdYOpt = NEWDARRAY (n * 2) ;

			for (nb = 0; nb < k; nb++)
			{
				fill (afin, p, 0) ;

				matmult (px, &n, &p, pr, &p, y) ;

				{			// OrderArray für die Y Matrix nach den SD der Spalten ordnen
					colsd (y, &n, &p, pdYColSd) ;
					order_decr (pdYColSd, &p, pnOrder) ;
					
					double *pdYOrdered = NEWDARRAY (n * p) ;
					SetColOrder (y, &n, &p, pnOrder, pdYOrdered) ;
		
					free (y) ;
					y = pdYOrdered ;
				}

				double dFoo ;
					//	res <- gridplane(y[,1:2],m,method=method)
					//	afin[1:2] <- res$alphamax
				gridplane (y, &n, &m, &method, &dFoo, afin) ;

					

									//	pdYOpt müsste eigentlich nur n x 1 sein, später wird jedoch ein temporäres colbind (+1 col)
									//	gemacht, deshalb wird hier schon eine col mehr allociert

					//	yopt <- y[,1:2]%*%res$alphamax

				int nTwoCols = 2 ;
				int nOneCol = 1 ;
				matmult (y, &n, &nTwoCols, afin, &nOneCol, pdYOpt) ;

				int j ;

				for (j = 2; j < p; j++)
				{
					double adAlphaMax[2] ;
						//	cbind(yopt,y[,j])
					matcpy (pdYOpt + n, y + j*n, n) ;		//	-> in pdYOpt + n steht y[,j] drinnen

						//	res <- gridplane(cbind(yopt,y[,j]),m,method=method)
					gridplane (pdYOpt, &n, &m, &method, 0, adAlphaMax) ;
//					double donefoo = 1 ;
//					gridplane2 (pdYOpt, &n, afin, &m, &donefoo, &method, 0, adAlphaMax) ;

						//	y[,j]*res$alphamax[2]
					VectorMult (pdYOpt + n, &n, adAlphaMax[1]) ;

						//	afin[1:(j-1)] <- afin[1:(j-1)]*res$alphamax[1]

						//	yopt <- yopt*res$alphamax[1]+y[,j]*res$alphamax[2]
							//	yopt <- yopt*res$alphamax[1]
					VectorMult (pdYOpt, &n, adAlphaMax[0]) ;

							//	yopt <- yopt+y[,j]*res$alphamax[2]
					VectorAddVector (pdYOpt, &n, pdYOpt + n) ;

					VectorMult (afin, &j, adAlphaMax[0]) ;

					//	afin[j] <- res$alphamax[2]
					afin[j] = adAlphaMax [1] ;
				}

			//	verified...

				double dObjf ;

				double *pd_y_x_afin = NEWDARRAY (n);
				matmult (y, &n, &p, afin, &nOneCol, pd_y_x_afin) ; 

				if (method)
					sd  (pd_y_x_afin, &n, &dObjf) ;
				else
					mad  (pd_y_x_afin, &n, &dObjf) ;

				double dObjF_old = dObjf, dObjF_best = dObjf ;

				matcpy (pdAfinBest, afin, p) ;

				int i ;
				for (i = 0; i < maxiter; i++)
				{
					for (j = 0; j < p; j++)
					{
						double adAlphaMax[2] ;

						matcpy (pdYOpt + n, y + n * j, n) ;
						double dDiv = pow (2, i+1) ;
						gridplane2 (pdYOpt, &n, afin + j, &m, &dDiv, &method, NULL, adAlphaMax) ;

						{	//	yopt <- yopt*res$alphamax[1]+y[,j]*res$alphamax[2]
								//	y[,j]*res$alphamax[2]
							VectorMult (pdYOpt + n, &n, adAlphaMax[1]) ;

								//	afin[1:(j-1)] <- afin[1:(j-1)]*res$alphamax[1]
								
									//	yopt <- yopt*res$alphamax[1]
							VectorMult (pdYOpt, &n, adAlphaMax[0]) ;

									//	yopt <- yopt+y[,j]*res$alphamax[2]
							VectorAddVector (pdYOpt, &n, pdYOpt + n) ;
						}

						VectorMult (afin, &p, adAlphaMax[0]) ;
						afin [j] += adAlphaMax[1] ;

					}
					double dAfinSqrtSumSqrt ;
					
					{	//	afin=afin/sqrt(sum(afin^2))
						dAfinSqrtSumSqrt= calcsqrtsumsqr (afin, p) ;
						
						VectorMult (afin, &p, 1 / dAfinSqrtSumSqrt) ;
					}

						//	y%*%afin = (n x p) X (p x 1) -> n x 1
					matmult (y, &n, &p, afin, &nOneCol, pd_y_x_afin) ; 

					if (method)
						sd  (pd_y_x_afin, &n, &dObjf) ;
					else
						mad  (pd_y_x_afin, &n, &dObjf) ;

					if (dObjf > dObjF_old)
					{
							dObjF_best = dObjf ;
							dObjF_old = dObjf ;

							matcpy (pdAfinBest, afin, p) ;
							VectorMult (pdAfinBest, &p, 1 / calcsqrtsumsqr (afin, p)) ;
					}
				}

				pObj [nb] = dObjF_best * dObjF_best ;
				SetVectorOrder (pdAfinBest, &p, pnOrder, l + nb * p) ;

				{	// pr <- pr-l[,nb]%*%t(l[,nb])  
					double *pd_lnb_x_t_lnb = NEWDARRAY (p * p ) ;		//	p x p
					matmult (l + nb * p, &p, &nOneCol, l + nb * p, &p, pd_lnb_x_t_lnb) ;
					MatSubMat (pr, pd_lnb_x_t_lnb, &p, &p) ;
				}

				free (pd_y_x_afin) ;

//				free (pdYOpt) ;

			}

			free (pdYColSd) ;
			free (pdAfinBest) ;

			free (pdYOpt) ;
			free (pnOrder) ;
			free (afin) ;
			free (y) ;
			free (pr) ;

			//	z <- x %*% l
			matmult (px, &n, &p, l, &k, z) ;
		}

//		free (z) ;
//		free (l) ;

	}
*/
/////////////////
//  PP - Code  //
/////////////////
/*
l1median <-
function(X)
{
	maxstep <- 200
	ittol <- 10^(-3)
	m <- apply(X, 2, median)
	k <- 1
	while(k < maxstep)
	{
		mold <- m
		Y <- scale(X, center = m, scale = F)
                #Ynorm <- apply(Y, 1, vecnorm)
                Ynorm <- sqrt(apply(Y^2, 1, sum))
		w <- ifelse(Ynorm == 0, 0, 1/Ynorm)
		W <- diag(w)
		delta <- apply(W %*% Y, 2, sum)/sum(w)
                #nd <- vecnorm(delta)
                nd <- sqrt(sum(delta^2))
		if(nd < ittol)
			maxhalf <- 0
		else maxhalf <- log(nd/ittol)/log(2)
		m <- mold
		nstep <- 0
		while((obj(X, m) >= obj(X, mold)) && (nstep <= maxhalf)
			) {
			nstep <- nstep + 1
			m <- mold + (delta/(2^nstep))
		}
		if(nstep > maxhalf)
			return(mold)
		k <- k + 1
	}
	if(k > maxstep)
		stop("iteration failed")
}
*/

/*
"obj" <-
function(X, m)
{
	//	Y skalieren
	//	jedes Element von Y ^2 & zeilenweiße die Summe bilden dann wurzek
	//	alle wurzeln aufsummieren


	Y <- scale(X, center = m, scale = F)
        #Ynorm <- apply(Y, 1, vecnorm)
    Ynorm <- sqrt(apply(Y^2, 1, sum))

	S <- sum(Ynorm)
	return(S)
}
*/

	extern void obj (double *pdData, int *pnRows, int *pnCols, double *pdm, double *pds)
	{		//	validated
		int i, j ;
		int nRows = *pnRows ;
		*pds = 0 ;
		double dCurSum ;
		double dCurValue ;
		for (i = nRows - 1; i >= 0; i--)
		{
			dCurSum = 0 ;

			for (j = *pnCols - 1; j >= 0; j--)
			{
				dCurValue = pdData [j * nRows + i] - pdm[j] ;
				dCurSum += dCurValue * dCurValue ;
			}
			*pds += ::sqrt (dCurSum) ;
		}
	}
	
	void GetSqrtofRowSumsofSquares (double *pdData, int *pnRow, int *pnCol, double *pdRet)
	{
		int i, j ;
		
		int nRow = *pnRow ;
		for (i = nRow - 1; i >= 0; i--)
		{
			pdRet [i] = 0 ;
			for (j = *pnCol - 1; j >= 0; j--)
			{
				double &dCurVal = pdData [j * nRow + i] ;
				pdRet [i] += dCurVal * dCurVal ;
			}
			pdRet [i] = ::sqrt ((double) (pdRet [i])) ;
		}
	}

/////////////////////////
// new l1median start  //
/////////////////////////


	inline double AddSqr (double &d, double pAdd)
	{
		return d += pAdd * pAdd ; 
	}

	class L1MinStruct
	{
	public:
		L1MinStruct (int n, int p, double *pdX, double *pdParScale = 0) ;
		~L1MinStruct () ;

		int m_n, m_p, m_pn ;
		double *m_pdX, *m_pdX_, *m_pdDi, *m_pdM, *m_pdParScale ; //, *m_pdM ;

		double calcall (double *pdM, double *pdMRet) ;
		double calObj (double *pdM) ;
		//void sqrtrowsumssq () ;

		int m_nCFn, m_nCGr ;

	} ;

	L1MinStruct::L1MinStruct (int n, int p, double *pdX, double *pdParScale) : m_n (n), m_p (p), m_pn (p * n), m_pdX (pdX), m_pdParScale (pdParScale)
	{
		m_pdX_ = new double [n * p] ;
		m_pdDi = new double [n] ;
		m_pdM = new double [p] ;	//	set m_pdM to some strange value!!
		m_nCFn = m_nCGr = 0 ;
	}

	L1MinStruct::~L1MinStruct ()
	{
		delete [] m_pdX_ ;
		delete [] m_pdDi ;
		delete [] m_pdM ;
	}

	void Rtprintf (double *pd, int n)
	{
		for (; n; n--)
			Rprintf ("%0.10f\t", *pd++) ;
		Rprintf ("\r\n") ;
	}

	double L1MinStruct::calcall (double *pdM, double *pdMRet)
	{
		m_nCGr++ ;
/*
		memcpy (m_pdM, pdM, m_p * sizeof (double)) ;

		if (m_pdParScale)
			VectorMultVector (m_pdM, &m_p, m_pdParScale) ;
*/
//Rprintf ("grad:\t") ;
//Rtprintf (pdM, m_p) ;
		int r, c ;
		double *pdX = m_pdX + m_pn, *pdX_ = m_pdX_ + m_pn ;

//		double *pdX = m_pdX, *pdX_ = m_pdX_ ;

		for (r = m_n - 1; r != -1; r--)
			m_pdDi[r] = 0 ;


        //    X. <- centr(X,m)
         //   d.i <- rowSums(X.^2)
		
		for (c = m_p - 1; c != -1; c--)
		{
			double dXm = m_pdParScale ? pdM[c] * m_pdParScale[c] : pdM[c] ;
			//for (r = 0; r < m_n; r++)
			for (r = m_n - 1; r != -1; r--)
			{
				*--pdX_ = *--pdX - dXm ;
				//*pdX_++ = *pdX++ - dXm ;
				m_pdDi[r] += *pdX_ * *pdX_ ;



			}
		}

		pdX_ = m_pdX_ + m_pn ;

		for (r = m_n - 1; r != -1; r--)
		{
			m_pdDi[r] = ::sqrt ((double) m_pdDi[r]) ;

		}

		//	MED <- colSums(X. / sqrt (d.i))


		for (c = m_p - 1; c != -1; c--)
		{
			double &dXm = pdMRet[c] = 0 ;
			for (r = m_n - 1; r != -1; r--)
				dXm -= *--pdX_ / m_pdDi[r] ;
		}

		return 0 ;
	}

	double L1MinStruct::calObj (double *pdM)
	{
		m_nCFn++ ;
//Rprintf ("obj: \t") ;
//Rtprintf (pdM, m_p) ;

		memcpy (m_pdM, pdM, m_p * sizeof (double)) ;

		if (m_pdParScale)
			VectorMultVector (m_pdM, &m_p, m_pdParScale) ;

		//double *pdCur = m_pdX +  m_n * m_p ;

		double dSum = 0, dRowSum ;
		int r, c ;
		for (r = m_n - 1; r != -1; r--)
		{
			dRowSum = 0 ;
			for (c = m_p - 1; c != -1; c--)
				AddSqr (dRowSum, m_pdX [r + c * m_n] - m_pdM [c]) ;

			dSum += ::sqrt (dRowSum) ;
		}
//Rprintf ("Obj. val: %f\r\n", dSum) ;
		return dSum ;
	}

	void l1objg(int n, double *pdCurCenter, double *pdRetCenter, void *pDat)
	{
		((L1MinStruct *) pDat)->calcall (pdCurCenter, pdRetCenter) ;
	}

	void l1obj_(int n, double *pdCurCenter, double *ddObj, void *pDat)
	{
		*ddObj = ((L1MinStruct *) pDat)->calObj (pdCurCenter) ;
	}

	extern void l1objg(int n, double *pdCurCenter, double *pdRetCenter, void *pDat) ;
	extern void l1obj_(int n, double *pdCurCenter, double *ddObj, void *pDat) ;

	double l1obj (int n, double *pdCurCenter, void *pDat)
	{
		return ((L1MinStruct *) pDat)->calObj (pdCurCenter) ;
	}

	extern void l1median_NM(int *pnParam, double *pdParam, double *pdData, double *pdParScale, double *pdMRet)
	{
		int &n = pnParam [0], &p = pnParam [1], &nMaxStep = pnParam[2], &nFail = pnParam [3], &nTrace = pnParam [4], &nFnCount = pnParam [5] ;
		double &dAbsTol = pdParam[0], &dRelTol = pdParam[1], &dRet = pdParam[2], &dAlpha = pdParam[3], &dBeta = pdParam[4], &dGamma = pdParam[5] ;

		double *pdStart = new double [p] ;

		memcpy (pdStart, pdMRet, sizeof (double) * p) ;

		VectorDivVector (pdStart, &p, pdParScale) ;

		L1MinStruct minstruc (n, p, pdData, pdParScale) ;

		nmmin (p, pdStart, pdMRet, &dRet, l1obj, &nFail, dAbsTol, dRelTol, (void *) &minstruc, dAlpha, dBeta, dGamma, nTrace, &nFnCount, nMaxStep) ;

		VectorMultVector (pdMRet, &p, pdParScale) ;

		delete [] pdStart ;
	}

	extern void l1median_CG(int *pnParam, int *pnParam_Out, double *pdParam, double *pdParam_Out, double *pdData, double *pdParScale, double *pdMRet)
	{
		int &n = pnParam [0], &p = pnParam [1], &nMaxStep = pnParam[2], &nTrace = pnParam [3], nType = pnParam [4] ;
		int &nFail = pnParam_Out [0], &nFnCount = pnParam_Out [1], &nGrCount = pnParam_Out [2] ;

		double &dAbsTol = pdParam[0], &dRelTol = pdParam[1] ;
		double &dRet = pdParam_Out [0] ;
		double *pdStart = new double [p] ;

		memcpy (pdStart, pdMRet, sizeof (double) * p) ;

		VectorDivVector (pdStart, &p, pdParScale) ;

		L1MinStruct minstruc (n, p, pdData, pdParScale) ;

		cgmin (p, pdStart, pdMRet, &dRet, l1obj, l1objg, &nFail, dAbsTol, dRelTol, (void *) &minstruc, nType, nTrace, &nFnCount, &nGrCount, nMaxStep) ;

		VectorMultVector (pdMRet, &p, pdParScale) ;

		delete [] pdStart ;
	}

	extern void l1median_BFGS (int *pnParam_In, int *pnParam_Out, double *pdParam_In, double *pdParam_Out, double *pdData, double *pdParScale, double *pdMRet)
	{
		int &n = pnParam_In [0], &p = pnParam_In [1], &nMaxStep = pnParam_In[2], &nTrace = pnParam_In [3], nReport = pnParam_In [4] ;
		int &nFail = pnParam_Out [0], &nFnCount = pnParam_Out [1], &nGrCount = pnParam_Out [2] ;

		double &dAbsTol = pdParam_In[0], &dRelTol = pdParam_In[1] ;
		double &dRet = pdParam_Out [0] ;

		VectorDivVector (pdMRet, &p, pdParScale) ;

		L1MinStruct minstruc (n, p, pdData, pdParScale) ;

		int *pnMask = new int [p] ;
		for (int i = p - 1; i != -1; i--)
			pnMask[i] = 1 ;

		vmmin (p, pdMRet, &dRet, l1obj, l1objg, nMaxStep, nTrace, pnMask, dAbsTol, dRelTol, nReport, (void *) &minstruc, &nFnCount, &nGrCount, &nFail) ;

		VectorMultVector (pdMRet, &p, pdParScale) ;

		delete [] pnMask ;
	}

	extern void l1median_SA (int *pnParam_In, int *pnParam_Out, double *pdParam_In, double *pdParam_Out, double *pdData, double *pdParScale, double *pdMRet)
	{
		int &n = pnParam_In [0], &p = pnParam_In [1], &nMaxStep = pnParam_In[2], &nTrace = pnParam_In [3], nTMax = pnParam_In [4] ;
		int &nFnCount = pnParam_Out [0] ;

		double &dTempInit = pdParam_In [0] ;
		double &dRet = pdParam_Out [0] ;

		VectorDivVector (pdMRet, &p, pdParScale) ;

		L1MinStruct minstruc (n, p, pdData, pdParScale) ;

		samin (p, pdMRet, &dRet, l1obj, nMaxStep, nTMax, dTempInit, nTrace, (void *) &minstruc) ;
		nFnCount = minstruc.m_nCFn ;

		VectorMultVector (pdMRet, &p, pdParScale) ;
	}
	
	extern void l1median_NLM (int *pnParam, double *pdParam, double *pdData, double *pdMRet, double *pdTypSize)
	{
		int &n = pnParam [0], &p = pnParam [1], &nMaxStep = pnParam[2], &nFail = pnParam [3], /*&nTrace = pnParam [4],*//* &nFnCount = pnParam [5],*/ &nMsg = pnParam[6] ;
		double &dTol = pdParam[0], &dRet = pdParam[1] ;

		double *pdStart = new double [p] ;

		memcpy (pdStart, pdMRet, sizeof (double) * p) ;

		L1MinStruct minstruc (n, p, pdData) ;

		nMsg = 8 ;

		double *pdRetGrad = new double [p] ;
			double *pdWrk = new double [p * 8] ;
		double *pdA = new double [p * p] ;
		optif9 (	p,
					p,
					pdStart,
					l1obj_,
					l1objg,
					NULL,
					(void *) &minstruc,
					pdTypSize,
					1,			/*	fscale	*/
					1,			/*	method	*/
					1,			/*	iexp	*/
					&nMsg,
					-1,			/*	ndigit	*/	//	// C: -1 / R: 12
					nMaxStep,
					1,			/*	iagflag	*/
					0,			/*	iahflag	*/
					1e-6,		/*	dlt		*/		// 1, //1e-6,
					dTol, 		/*	gradtl	*/		//pow(DBL_EPSILON, 0.25), //0.1,
					1000,		/*	stepmx	*/		//	C:1e-6 / R: 1000
					dTol,		/*	steptl	*/		//sqrt(DBL_EPSILON), //1e-6,			
					pdMRet,	
					&dRet,
					pdRetGrad,
					&nFail,
					pdA,
					pdWrk,
					&nMaxStep) ;

		delete [] pdStart ;
		delete [] pdRetGrad ;
		delete [] pdWrk ;
		delete [] pdA ;
	}

///////////////////////
// new l1median end  //
///////////////////////

/////////////////////////
// old l1median start  //
/////////////////////////

	extern void l1median(double *pdData, int *pnRow, int *pnCol, double *pdMRet, int *pnRet, int *pnMaxStep, double *pdItTol)
	{		//	abgesehen von (k > nMaxStep) (am schluss) validated

		*pnRet = 1 ;
		int nMaxStep  = pnMaxStep ? *pnMaxStep: 200 ;
		double dItTol = pdItTol ? *pdItTol: 0.00000001 ;

			//	m <- apply(X, 2, median)
		int nDataCol = *pnCol ;
		int nDataRow = *pnRow ;
		int i,j ;
		double *pdM = NEWDARRAY (nDataCol) ;
	//	double *pdMold = NEWDARRAY (nDataCol) ;
		int nDataSize = nDataRow * nDataCol ;
		double *pdY = NEWDARRAY (nDataRow * nDataCol) ;		// nRow x nCol
		double *pdW = NEWDARRAY (nDataRow * nDataRow) ;
		double dSumW ;										//	nRow x nRow

		double *pdWxY = NEWDARRAY (nDataRow * nDataCol) ;
		double *pdDelta = NEWDARRAY (nDataCol) ;
		double dND ;
		double dMaxHalf = 0 ;

		ColMedian (pdData, pnRow, pnCol, pdM, TRUE) ;
		int k = 1 ;

		while(k < nMaxStep)
		{
			matcpy (pdMRet, pdM, nDataCol) ;

				//Y <- scale(X, center = m, scale = F)
			matcpy (pdY, pdData, nDataSize) ;
			MatSubstractVectorByRow (pdY, pnRow, pnCol, pdM) ;



				//#Ynorm <- apply(Y, 1, vecnorm)
			{		//Ynorm <- sqrt(apply(Y^2, 1, sum))
					//w <- ifelse(Ynorm == 0, 0, 1/Ynorm)
					//W <- diag(w)

				dSumW = 0 ;
				fill (pdW, nDataRow * nDataRow, 0) ;
				double dCurValue ;
				for (i = nDataRow - 1; i >= 0; i--)				//	ab hier in calcsquared rowsums in eigene fkt stecken (auch für obj?!)
				{
					double &dCurSum = pdW[i + i * nDataRow] ;
					dCurSum  = 0 ;

					for (j = nDataCol - 1; j >= 0; j--)
					{
						dCurValue = pdY [i + j * nDataRow] ;
						dCurSum += dCurValue * dCurValue ;
					}
					dCurSum = ::sqrt ((double) dCurSum) ;
					if (dCurSum != 0.0)
						dCurSum = 1 / dCurSum ;
					dSumW += dCurSum ;
				}
			}
				
			{		//delta <- apply(W %*% Y, 2, sum)/sum(w)
				matmult (pdW, pnRow, pnRow, pdY, pnCol, pdWxY) ;
						
				ColSums (pdWxY, pnRow, pnCol, pdDelta) ;
				VectorMult (pdDelta, pnCol, 1 / dSumW) ;
			}


			{       //	#nd <- vecnorm(delta)
					//	nd <- sqrt(sum(delta^2))
				dND = 0 ;
				for (i = nDataCol - 1; i >= 0; i--)
					dND += pdDelta [i] * pdDelta [i] ;
				dND = ::sqrt((double) dND) ;
			}

			if(dND < dItTol)
				dMaxHalf = 0 ;
			else
				dMaxHalf = ::log ((double) (dND/dItTol))/::log((double) 2) ;

					//	m <- mold
			matcpy (pdM, pdMRet, nDataCol) ;
			
			int nStep = 0 ;

			double dObjM ;
			double dObjMRet ;
			obj(pdData, pnRow, pnCol, pdM, &dObjM) ;
			obj(pdData, pnRow, pnCol, pdMRet, &dObjMRet) ;

			while(dObjM >= dObjMRet && (nStep <= dMaxHalf))
			{
				nStep ++ ;
					//	m <- mold + (delta/(2^nstep))
				double dPowStep = ::pow ((double) 2.0, nStep) ;
				for (i = nDataCol - 1; i >= 0; i--)
					pdM[i] = pdMRet[i] + pdDelta[i] / dPowStep ;

				obj(pdData, pnRow, pnCol, pdM, &dObjM) ;
				obj(pdData, pnRow, pnCol, pdMRet, &dObjMRet) ;
			}

			if(nStep > dMaxHalf)
				goto free_mem ;	//	return!!
			k ++ ;
		}

		if(k > nMaxStep)
			*pnRet = 0 ;

		free_mem:
		free (pdM) ;
	//	free (pdMold) ;
		free (pdY) ;
		free (pdW) ;
		free (pdWxY) ;
		free (pdDelta) ;
		return ;
	}

///////////////////////
// old l1median end  //
///////////////////////

/*

  prcomp.rob <-
function(X, k = 0, sca = "A", scores = F)
{
#
# C. Croux (Brussels, Belgium), P. Filzmoser (Vienna, Austria)
#
# k is the number of Principal Components that will be computed
# sca is the scale estimator that will be use as projection index:
#   "mad", "tau", and "A" estimators can be chosen
# 
# Output:
#     $scale .... vector of first k eigenvalues
#     $rotation ... matrix with the first k eigenvectors in the columns
#     $scores ... matrix with the first k score vectors in the columns (if scores=T)
# 
	n <- nrow(X)
	p <- ncol(X)
	if(k == 0)
		p1 <- min(n, p)
	else
		p1 <- k
	S <- rep(1, p1)
	V <- matrix(1:(p * p1), ncol = p1, nrow = p)
	P <- diag(p)	
m <- l1median(X)	
#	m <- apply(X, 2, median)
	Xcentr <- scale(X, center = m, scale = F)
	for(k in 1:p1) {
		B <- Xcentr %*% P
                #Bnorm <- apply(B, 1, vecnorm)
                Bnorm <- sqrt(apply(B^2,1,sum))
		A <- diag(1/Bnorm) %*% B
		Y <- A %*% P %*% t(X)
		if(sca == "mad")
			s <- apply(Y, 1, mad)
		if(sca == "tau")
			s <- apply(Y, 1, scale.tau)
		if(sca == "A")
			s <- apply(Y, 1, scale.a)
		j <- order(s)[n]
		S[k] <- s[j]
		V[, k] <- A[j,  ]
		if(V[1, k] < 0)
			V[, k] <- (-1) * V[, k]
		P <- P - (V[, k] %*% t(V[, k]))
	}
	if(scores)
	{
		list(scale = S, rotation = V, scores = Xcentr %*% V)
	}
	else
		list(scale = S, rotation = V)
}
*/


/*	extern void prcomp_rob (double *pdData, int *pnParams, double *pdS, double *pdV, double *pdScores)
	{
//		REprintf ("Check this out") ;
		int &n = pnParams[0] ;
		int &p = pnParams[1] ;
		int &k = pnParams[2] ;
//		int &nScale = pnParams[3] ;
		int &nScores = pnParams[4] ;
		int p1 ;
		if(k)
			p1 = k;
		else
			p1 = n < p ? n : p ;// min(n, p) ;
		//		V <- matrix(1:(p * p1), ncol = p1, nrow = p)
		//	pdV		p x p1
		int i ;
		for (i = p * p1 - 1; i >= 0; i--)
			pdV[i] = i + 1 ;

		//	P <- diag(p)
		double *pdP = new double [p * p] ;			//	P: p x p
		fillDiag (pdP, &p, &p) ;
		
		double *pdM = new double [p] ;
		int nRet ;
		//	m <- l1median(X)	
		l1median (pdData, &n, &p, pdM, &nRet, 0, 0) ;
		
		double *pdXcentr = new double [n * p] ;		//	pdXcentr: n x p
		matcpy (pdXcentr, pdData, n * p) ;
		//	Xcentr <- scale(X, center = m, scale = F)
		MatSubstractVectorByRow (pdXcentr, &n, &p, pdM) ;
		
		double *pdB = new double [n * p] ;			//	pdB: n x p
		double *pdBnorm = new double [n] ;
		double *pdA = new double [n * p] ;			//	pdA: n x p

		double *pdtData = new double [p * n] ;		//	pdtData p x n
		t(pdData, &n, &p, pdtData) ;

		double *pdTemp = new double [p * n] ;		//	for 	Y <- A %*% P %*% t(X)  (n x p) * (p x p) * (p x n)
		double *pdY = new double [n * n] ;			//	pdY: n x n
		double *pds = new double [n] ;				//	pds: n

		double *pdV2 = new double [p * p] ;			//	pdV2: p x p		enthält pdV²
		int	*pnOrder = new int [n] ;
		int j ;

		for (i = 0; i < p1; i++)
		{
			//	S <- rep(1, p1)
			pdS[i] = 1 ;
			//	B <- Xcentr %*% P
			matmult (pdXcentr, &n, &p, pdP, &p, pdB) ;

			//Bnorm <- sqrt(apply(B^2,1,sum))
			GetSqrtofRowSumsofSquares (pdB, &n, &p, pdBnorm) ;
			InvertVector (pdBnorm, &n) ;

			//	A <- diag(1/Bnorm) %*% B		//	n x n %*% n x p  -> n x p
			matcpy (pdA, pdB, n * p) ;
			MatMultVectorByCol (pdA, &n, &p, pdBnorm) ;

			//	Y <- A %*% P %*% t(X)			//	n x p %*% p x p %*% p x n -> n x n
			matmult(pdA, &n, &p, pdP, &p, pdTemp) ;
			matmult (pdTemp, &n, &p, pdtData, &n, pdY) ;

//			if (!nScale)
				//	s <- apply(Y, 1, mad)
				rowmad (pdY, &n, &n, pds) ;
//			else if (nScale == 1)
//				s <- apply(Y, 1, scale.tau)
//			else
//				s <- apply(Y, 1, scale.a)
			order (pds, &n, pnOrder) ;
			j = pnOrder [n - 1] ;
			pdS[i] = pds[j] ;
			//	V[, k] <- A[j,  ]
			GetRow (pdA, &n, &p, &j, pdV + p * i) ;

			if(pdV[p * i] < 0)
//				V[, k] <- (-1) * V[, k]
				VectorMult (pdV + p * i, &p, -1) ;

//			P <- P - (V[, k] %*% t(V[, k]))		p x 1 x 1 x p
			j = 1 ;
			matmult (pdV + p * i, &p, &j, pdV + p * i, &p, pdV2) ;
			j = p * p ;
			VectorSubstractVector (pdP, &j, pdV2) ;

		}

		if (nScores)
			matmult (pdXcentr, &n, &p, pdV, &p1,pdScores) ;		//	n x p * p x p1

		delete [] pdP ;
		delete [] pdM ;
		delete [] pdB ;
		delete [] pdBnorm ;
		delete [] pdtData ;

		delete [] pdTemp ;
		delete [] pdY ;
		delete [] pds ;
		delete [] pnOrder ;
		delete [] pdV2 ;
			

	}
*/
