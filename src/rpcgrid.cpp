#include "math.h"
#include "R.h"
#include "rsubst.h"
#include "fastpca.h"



extern 	void rpcgrid (double *px, int *pnParams, double *pObj, double *pLoadings, double *pScores)
	{
		int n  = pnParams [0] ;		//	nRows
		int p  = pnParams [1] ;		//	nCols

		int k	= pnParams [2] ;
		int m	= pnParams [3] ;
		int maxiter = pnParams [4] ;
		int	method = pnParams [5] ;
		int nAngleHalving = pnParams [6] ;
		int nScores = pnParams [7] ;
		int n2DimFact = pnParams [8] ;

		dev_Type *pfDev = GetDevFunction (method) ;

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

//		double *z = pScores /*= (double *) malloc (n * k * sizeof (double))*/ ;		//	n x k	//	matrix with scores
		if (nScores)		
			fill (pScores, n * k, R_NaReal) ;

		double *l = pLoadings /*= (double *) malloc (p * k * sizeof (double))*/ ;		//	p x k	//	matrix with loadings

		fill (pLoadings, p * k, R_NaReal) ;

		if (p == 2)
		{
			int nCurM = m * n2DimFact ;
			gridplane (px, &n, &nCurM, pfDev, pObj, l) ;

			if (k > 1)
			{
				l[0 + 1 * p] = -l[1] ;
				l[1 + 1 * p] = l[0] ;

				double *pFatMat = NEWDARRAY (n) ;		//	n x 1

				int nOneCol = 1 ;
				matmult (px, &n, &p, l + p, &nOneCol, pFatMat) ;

/*				if (method)
					colsd (pFatMat, &n, &nOneCol, pObj + 1) ;
				else
					colmad (pFatMat, &n, &nOneCol, pObj + 1) ;*/
				coldev (pFatMat, &n, &nOneCol, pObj + 1, pfDev) ;

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
				gridplane (y, &n, &m, pfDev, &dFoo, afin) ;



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
						//        %%% this is the difference between rpcgrid.R & rpcgrid2.R
						//	res <- gridplane2(cbind(yopt,y[,j]),afin[j],m,div=2^i,method=method)

//					double dDiv = pow (2, i) ;
//					gridplane2 (pdYOpt, &n, afin+ j, &m, &dDiv, &method, 0, adAlphaMax) ;
					gridplane (pdYOpt, &n, &m, pfDev, 0, adAlphaMax) ;

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

//				if (method)
//					sd  (pd_y_x_afin, &n, &dObjf) ;
//				else
//					mad  (pd_y_x_afin, &n, &dObjf) ;
				pfDev (pd_y_x_afin, &n, &dObjf) ;

				double dObjF_old = dObjf, dObjF_best = dObjf ;

				matcpy (pdAfinBest, afin, p) ;

				int i ;
				for (i = 0; i < maxiter; i++)
				{
					for (j = 0; j < p; j++)
					{
						double adAlphaMax[2] ;

						matcpy (pdYOpt + n, y + n * j, n) ;
						double dDiv = nAngleHalving ? 2<< i : 1 ;
						gridplane2 (pdYOpt, &n, afin + j, &m, &dDiv, pfDev, NULL, adAlphaMax) ;

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

//					if (method)
//						sd  (pd_y_x_afin, &n, &dObjf) ;
//					else
//						mad  (pd_y_x_afin, &n, &dObjf) ;
					pfDev (pd_y_x_afin, &n, &dObjf) ;

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

			//	pScores <- x %*% l
			if (nScores)
				matmult (px, &n, &p, l, &k, pScores) ;
		}

		VectorSqrt (pObj, &k, pObj) ;
	}
