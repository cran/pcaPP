#include "math.h"
#include "R.h"
#include "rsubst.h"
#include "fastpca.h"



extern 	void rpcgrid (double *pdx, int *pnParams, double *pdObj, double *pdl, double *pdScores)
{
		int n  = pnParams [0] ;		//	nRows
		int p  = pnParams [1] ;		//	nCols

		int k	= pnParams [2] ;
		int m	= pnParams [3] ;
		int maxiter = pnParams [4] ;
		int	method = pnParams [5] ;
		int nAngleHalving = pnParams [6] ;
		int nScores = pnParams [7] ;
		int n2DimFact = pnParams [8] * m;

		dev_Type *pfDev = GetDevFunction (method) ;

		// grid search for PCA, using Householder transformation for orthogonality of loadings
		//
		// x ... centered (!) data matrix
		// k ... number of components to determine
		// m ... number of directions to search
		// maxiter ... maximum number of iterations
		// method ... how to estimate the standard deviation
		//        should be: "sd", "mad"
		// fact2dim ... multiplication factor for m for 2-dimensional data
		//

		//	z <- matrix(NA,nrow=n,ncol=k) # matrix with scores
//		double *pdz = new double [n * k] ; // matrix with scores

		//	l <- matrix(0,nrow=p,ncol=k) # matrix with loadings
//		double *pdl = new double [p * k] ; // matrix with scores
		fill (pdl, p * k, 0) ;

		//	lorigin <- matrix(NA,nrow=p,ncol=k) # matrix with loadings in original space

		double *pdlOrigin = new double [p * k] ; // matrix with scores
		fill (pdlOrigin, p * k, R_NaReal) ;

		//	obj <- rep(NA,k) # vector with objective function
		//double *pdObj = new double [k] ; // matrix with scores
		//fill (pdObj, p * k, R_NaReal) ;

		int *pnOrder = new int [p] ; 

//		int nCurM = m * n2DimFact ;

		double *pdTempN = new double [n + p] ;
		double *pdTempP = pdTempN + n ;

		static int nTwo = 2 ;
		static int nOne = 1 ;

		int nFoo ;
		double dFoo ;

		

		if (p == 2)
		{
			// res <- gridplane(x,fact2dim*m,method=method)

			gridplane (pdx, &n, &n2DimFact, pfDev, pdObj, pdl) ;

			if (k > 1)
			{
				pdl[p] = -pdl[1] ;
				pdl[1 + p] = pdl[0] ;

				double *pFatMat = new double [n] ;		//	n x 1

				int nOneCol = 1 ;
				matmult (pdx, &n, &p, pdl + p, &nOneCol, pFatMat) ;

				coldev (pFatMat, &n, &nOneCol, pdObj + 1, pfDev) ;

				delete [] pFatMat ;
			}
		}
		else if (p > 2)
		{
						//	Transfo <- diag(p) # Transformation matrix
			double *pdTransfo = new double [p * p] ;
			double *pdTransfo_x_Edge = new double [p * p] ;
			fillDiag (pdTransfo, &p, &p) ;
						//	y <- x
			double *pdy = new double [n * p] ;
			double *pdY_Alloc = pdy ;
			double *pdU = new double [p * p] ;			// p x p

			matcpy (pdy, pdx, n * p) ;

			double *pdY_x_U = new double [n * p] ;		//	p x n
			double *pdY_x_U_Alloc = pdY_x_U ;


			

			double *pdAfinbest = new double [p] ;
			double pdAlphaMax [2] ;
			double dObjfbest ;//= new double [p * k] ;
			double *pdAfin = new double [p] ;			// p x 2
			double *pdYColSd = new double [p] ;			// p x 1
			double *pdYOrdered = new double [n * p] ;	// n x p
			double *pdYOpt	 = new double [n * 2] ;		// n x p

			double *pdBase = new double [p * p] ;		// p x p

			double *pdEdge = new double [p * p] ;		// p x p

			double dObjf, dObjfold ;
			double dnDiff ;
			int nCurDimCount ;								//	current number of dimensions

			double *pCurLPos ;

			int nb ;
			for (nb = 0; nb < k; nb++)
			{	  // loop over number of comp

				nCurDimCount = p - nb ;

				int p1 = p - nb ;//+ 1 ; // dimension will be reduced for subsequent PCs

				if (p1 == 1)
					fill (pdl + nb + nb * p, p - nb, 1) ;
				else if (p1 == 2)
				{
					gridplane (pdy, &n, &n2DimFact, pfDev, &dObjfbest, pdAfinbest) ;
					pdObj[nb] = dObjfbest ;
					pnOrder[0] = 0 ;
					pnOrder[1] = 1 ;
				}
				else
				{	// more than 2 dimensions from data left
					fill (pdAfin, p1, 0) ;
					{			// OrderArray für die Y Matrix nach den SD der Spalten ordnen
						coldev (pdy, &n, &nCurDimCount,pdYColSd,  pfDev) ;
						//colsd (pdy, &n, &p, pdYColSd) ;
						order_decr (pdYColSd, &nCurDimCount, pnOrder) ;
						SetColOrder (pdy, &n, &nCurDimCount, pnOrder, pdYOrdered) ;
					}


					gridplane (pdYOrdered, &n, &m, pfDev, &dFoo, pdAfin) ;
								//yopt <- yord[,1:2]%*%res$alphamax
					matmult (pdYOrdered, &n,  &nTwo, pdAfin, &nOne, pdYOpt) ;

					int j, i ;
					for (j = 2; j < p1; j++)
					{
						matcpy (pdYOpt + n, pdYOrdered + n * j, n) ;
						gridplane (pdYOpt, &n, &m, pfDev, &dObjf, pdAlphaMax) ;
						matmult (pdYOpt, &n, &nTwo, pdAlphaMax, &nOne, pdTempN) ;
						matcpy (pdYOpt, pdTempN, n) ;
						nFoo = j  ;
						VectorMult (pdAfin, &j, *pdAlphaMax) ;
						pdAfin[j] = pdAlphaMax [1] ;

					}

					dObjfbest = dObjfold = dObjf ;
					matcpy (pdAfinbest, pdAfin, 2) ;

					for (i = 0; i < maxiter; i++)
					{
						for (j = 0; j < p1; j++)
						{
							matcpy (pdYOpt + n, pdYOrdered + n * j, n) ;
							double dDiv = nAngleHalving ? 2<< i : 1 ;
							gridplane2 (pdYOpt, &n, pdAfin + j, &m, &dDiv, pfDev, &dObjf, pdAlphaMax) ;
							matmult (pdYOpt, &n, &nTwo, pdAlphaMax, &nOne, pdTempN) ;
							matcpy (pdYOpt, pdTempN, n) ;
							VectorMult (pdAfin, &p, *pdAlphaMax) ;
							pdAfin[j] += pdAlphaMax [1] ;
						}

									//	afin=afin/sqrt(sum(afin^2))
						VectorMult (pdAfin, &p, 1 / calcsqrtsumsqr (pdAfin, nCurDimCount)) ;

						if (dObjf > dObjfold)
						{
							dObjfold = dObjfbest = dObjf ;
										//	afinbest <- afin/sqrt(sum(afin^2))
							matcpy (pdAfinbest, pdAfin, p) ;
							VectorMult (pdAfinbest, &p, 1 / calcsqrtsumsqr (pdAfin, nCurDimCount)) ;
						} 
					}
					pdObj[nb] = dObjfbest ;

				}	// how many dimensions from data are left

				if (p1 == 1)
				{	// only 1 dimension of data left

								//	lorigin[,nb]=Transfo%*%l[,nb]
					matmult (pdTransfo, &p, &p, pdl + nb * p, &nOne, pdlOrigin + p * nb) ;
								//	if (method=="sd"){objflast <- sd(x%*%lorigin[,nb])} 
								//	else {objflast <- mad(x%*%lorigin[,nb])} 
								//	obj[nb] <- objflast

					matmult (pdx, &n, &p, pdlOrigin + p * nb, &nOne, pdTempN) ;
					pfDev (pdTempN, &n, pdObj + nb) ;
				}
				else
				{	// more than 1 dimension in data left
					pCurLPos = pdl + nb + nb * p ;
								// l[nord+nb-1,nb] <- afinbest # column nb of loadings matrix # !!!!!!!!!!!!!

					SetVectorOrder (pdAfinbest, &nCurDimCount, pnOrder, pCurLPos) ;

					fillDiag  (pdBase , &nCurDimCount, &nCurDimCount) ;
					dnDiff = calcsumsqrVectorDiff (pdBase, pCurLPos, p-nb) ;
					if (dnDiff > 1e-12)
					{
									//	if ((t(l[nb:p,nb])%*%Base[,1]) < 0)
						matmult (pCurLPos, &nOne, &nCurDimCount, pdBase, &nOne, &dFoo) ;
						if (dFoo < 0)
							VectorMult (pCurLPos, &nCurDimCount, -1) ;
									//	u=(Base[,1]-l[nb:p,nb])/sqrt(sum((Base[,1]-l[nb:p,nb])^2))
						
						matcpy (pdTempN, pdBase, nCurDimCount) ;
						VectorSubstractVector (pdTempN, &nCurDimCount, pCurLPos) ;			//	pdTempN = u  = 1 x nCurDimCount
						dFoo = ::sqrt (calcsumsqrVectorDiff (pdBase, pCurLPos, nCurDimCount)) ;
						VectorMult (pdTempN, &nCurDimCount, 1 / ::sqrt (calcsumsqrVectorDiff (pdBase, pCurLPos, nCurDimCount))) ;

									//	hlp=t(u)%*%Base
						matmult (pdTempN, &nOne, &nCurDimCount, pdBase, &nCurDimCount, pdTempP) ;	//	pdTempP = hlp = 1 x nCurDimCount

									//	U=Base - (t((hlp*2)) %*% u)
						VectorMult (pdTempP, &nCurDimCount, -2) ;
						matmult (pdTempP, &nCurDimCount, &nOne, pdTempN, &nCurDimCount, pdU) ;
						nFoo = nCurDimCount * nCurDimCount ;
						VectorAddVector (pdU, &nFoo, pdBase) ;
					}
					else
						matcpy (pdU, pdBase, nCurDimCount * nCurDimCount) ;

					// Back transformation to original pxp space:
///	Verified...
					if (!nb)
					{
						matcpy (pdlOrigin + p * nb, pdl + p * nb, p) ;
						matcpy (pdTransfo, pdU, p * p) ;
					}
					else
					{
						fillDiag (pdEdge, &p, &p) ;
						MatSetMat (pdEdge, pdU, p, nCurDimCount, nCurDimCount, nb, nb) ;
						matmult (pdTransfo, &p, &p, pdl + nb * p, &nOne, pdlOrigin + nb * p) ;

						matmult (pdTransfo, &p, &p, pdEdge, &p, pdTransfo_x_Edge) ;
						SwapPtrs (&pdTransfo_x_Edge, &pdTransfo) ;		//	memorybereiche von pdTransfo_x_Edge % pdTransfo swapen - damit sparen wir uns 1x kopieren

					}
					matmult (pdy, &n, &nCurDimCount, pdU, &nCurDimCount, pdY_x_U) ;
					SwapPtrs (&pdY_x_U, &pdy) ;
					pdy += n ;
				}
			}

			delete [] pdTransfo ;
			delete [] pdY_Alloc ;
			delete [] pdU ;
			delete [] pdY_x_U_Alloc ;

			delete [] pdAfin ;
			delete [] pdYColSd ;
			delete [] pdYOrdered ;
			delete [] pdYOpt ;
			delete [] pdBase ;
			delete [] pdAfinbest ;
			delete [] pdEdge ;
			delete [] pdTransfo_x_Edge ;

			order_decr (pdObj, &k, pnOrder) ;
			SetVectorOrder (pdObj, &k, pnOrder, pdTempP) ;
			matcpy (pdObj, pdTempP,k) ;
			SetColOrder (pdlOrigin, &p, &k, pnOrder, pdl) ;

		}


		if (nScores)
			matmult (pdx, &n, &p, pdl, &k, pdScores) ;
		
		delete [] pdlOrigin ;
		delete [] pdTempN ;
		delete [] pnOrder ;

}
