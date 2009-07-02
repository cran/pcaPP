#include "L1Median.h"
#include "ext.h"
#include "perftimer.h"

	inline double AddSqr (double &d, const double &pAdd)	{ return d += pAdd * pAdd ;  }

	double calObj (double *pdData, double *pdM, int n, int p)
	{
		double dSum = 0, dRowSum ;
		int r, c ;

		for (r = n - 1; r != -1; r--)
		{
			dRowSum = 0 ;

			for (c = p - 1; c != -1; c--)
				AddSqr (dRowSum, pdData [r + c * n] - pdM [c]) ;

			dSum += ::sqrt (dRowSum) ;
		}
		return dSum ;
	}

	void l1median_HoCr (int *pnParam_In, int *pnParam_Out, double *pdParam_In, double *pdParam_Out, double *pdData, double *pdMRet)
	{
		const double dLog2 = log (2) ;

		int &n = pnParam_In [0], &p = pnParam_In [1], &dwMaxit = pnParam_In[2], &dwTrace = pnParam_In [3] ;

		int &dwIterCount = pnParam_Out [0], &dwStepHalvingSteps = pnParam_Out [1] = 0, &dwCode = pnParam_Out[2] = 0, &nTime = pnParam_Out[3] ;
		double &dTol = pdParam_In[0] ;

		double &dObj = pdParam_Out [0], dObjOld, dND  ;
		IVecD vdMedian (p, pdMRet), vdMedianOld (p), vdNorms (n), vdNormsOrdered, vdWeights (n), vdDelta (p) ;
		IMatD mdData (n, p, pdData), mdDataZ (n, p) ;


		IVecDW vinter (n) ;

		int &k = dwIterCount, i, dwNStep, dwMaxHalf ;

		int nNHalf = n >> 1 ;

		CPerfTimer tim ;

		for (k = 0; k < dwMaxit; k++)
		{
			vdMedianOld << vdMedian ;

			dObjOld = k ? dObj : calObj (pdData, pdMRet, n, p) ;

			FC_ElOp<FC::FC_minus, double>::OpMV_row (mdData, vdMedian, mdDataZ) ;	//	mdDataZ =  scale (mdData, center = vdMedian, scale = FALSE)	

			SqrtrowSumSqs (mdDataZ, vdNorms) ;

			vdNorms.order (vinter) ;

			vdNormsOrdered = vdNorms (ISub (vinter)) ;

			mdData = mdData (ISub (vinter), ISub ()) ;
			mdDataZ = mdDataZ (ISub (vinter), ISub ()) ;

			double dSumWeights = 0 ;
			DWORD dwCountZeroDiff = 0 ;
			for (i = n - 1; i != (int) -1; i--)
			{
				double &dCurVal = vdNormsOrdered (i) ;

				if (dCurVal == 0)
				{
					vdWeights (i) = 0 ;
					dwCountZeroDiff++ ;

					if (i > nNHalf)	//	there's more than half of the values concentrated at the current median estimation -> return this estimation
					{
						if (dwTrace >= 1)
							Rprintf ("A concentration of %d >= n/2 = %d observations in one point has been detected (return code 3).\r\n", i + 1, nNHalf) ;
						dwCode = 3 ;
						dObj = dObjOld ;
						nTime = tim.GetTimeMS () ;
						return ;
					}


				}
				else
					dSumWeights += (vdWeights (i) = 1 / dCurVal) ;
			}

			ColSumWeighted (mdDataZ, vdWeights, vdDelta) ;
			vdDelta /= dSumWeights ;
			sqrtsumsq (vdDelta, dND) ;

			if (dwTrace >= 3)	Rprintf ("nd at %g in iteration %d (tol at %g)\r\n", dND, k, dTol) ;

			vdMedian += vdDelta ;

			if (dND < dTol)	//	converged
			{
				dObj = calObj (pdData, pdMRet, n, p) ;	//	calc current objective function
				nTime = tim.GetTimeMS () ;
				return  ;
			}

			dwMaxHalf = (int) ceil (log (dND / dTol) / dLog2) ;

			for (dwNStep = 0; dwNStep < dwMaxHalf && (dObj = calObj (pdData, pdMRet, n, p)) >= dObjOld; dwNStep++)
			{
				dwStepHalvingSteps++ ;
				vdDelta /= 2.0 ;

				FC_ElOp<FC::FC_plus, double>::OpVV  (vdMedianOld, vdDelta, vdMedian) ;
			}

			if (dwNStep >= dwMaxHalf)
			{	
				if (dwTrace >= 1)
					Rprintf ("step halving failed in %d steps (return code 2)\r\n", dwMaxHalf) ;
				dwCode = 2 ;
				vdMedian << vdMedianOld ;
				dObj = dObjOld ;
				nTime = tim.GetTimeMS () ;
				return ;
			}

		}
		nTime = tim.GetTimeMS () ;

		dwCode = 1 ;	//	algorithm did not converge

		if (dwTrace >= 1)
			Rprintf ("Algorithm did not converve (return code 1).\r\n") ;
	}

	void foo ()
	{
		IDim mydim(20) ;

		mydim.GetSize () ;
	}
