#include "math.h"
#include <stdlib.h>
#include "R.h"
#include <memory.h>
#include "fastpca.h"

#include "rsubst.h"


	void VectorAdd (double *x, int *pnSize, double dAdd)
	{
		int i ;
		for (i = *pnSize - 1; i >=0; i--)
			x[i] += dAdd ;
	}

	void VectorAddVector (double *x, int *pn, double *pdAdd)
	{
		int i ;
		for (i = *pn - 1; i >=0; i--)
			x[i] += pdAdd[i] ;
	}

	void VectorSubstractVector (double *x, int *pn, double *pdSubstract)
	{
		int i ;
		for (i = *pn - 1; i >=0; i--)
			x[i] -= pdSubstract[i] ;
	}

	void VectorMultVector (double *x, int *pn, double *pdMult)
	{
		int i ;
		for (i = *pn - 1; i >=0; i--)
			x[i] *= pdMult[i] ;
	}

	void MatDivVectorByColIdx (double *pdData, int n, int p, double *pdDiv, int *pdIdx)
	{
		int i ;
		for(i = p; i > 0; i--)
		{
			VectorDivVectorIdx (pdData, &n, pdDiv, pdIdx) ;
			pdData += n ;
		}
	}

	void VectorDivVectorIdx (double *x, int *pn, double *pdDiv, int *pnIdx)
	{
		int i ;
		for (i = *pn - 1; i >=0; i--)
			if (pnIdx[i])
				x[i] /= pdDiv[i] ;
	}

	void VectorDivVector (double *x, int *pn, double *pdDiv)
	{
		int i ;
		for (i = *pn - 1; i >=0; i--)
			x[i] /= pdDiv[i] ;
	}

	void MatAddMat (double *pMat1, double *pMat2, int *pnRow, int *pnCol)
	{
		int i; 
		for (i = *pnRow * *pnCol - 1; i >= 0; i--)
			*(pMat1++) += *(pMat2++) ;
	}

	void MatAddVectorByRow (double *pMat1, int *pnRow, int *pnCol, double *pAdd)
	{
		int i,j ; 

		for (j = *pnCol - 1; j >= 0; j--)
		{
			for (i = *pnRow - 1; i >= 0; i--)
				*(pMat1++) += *pAdd ;
			pAdd++ ;
		}
	}

	void MatMultVectorByCol (double *pdData, int *pnRow, int *pnCol, double *pMult)
	{
		int i ; 
		pdData += *pnCol * *pnRow  - 1;
		for (i = *pnCol * *pnRow - 1; i >= 0; i--)
			*(pdData--) *= pMult[i % *pnRow] ;
	}

	void MatDivVectorByCol (double *pdData, int *pnRow, int *pnCol, double *pDiv)
	{
		int i ; 
		pdData += *pnCol * *pnRow  - 1;
		for (i = *pnCol * *pnRow - 1; i >= 0; i--)
			*(pdData--) /= pDiv[i % *pnRow] ;
	}

	void MatSubstractVectorByRow (double *pMat1, int *pnRow, int *pnCol, double *pSubstract)
	{
		int i,j ; 

		for (j = *pnCol - 1; j >= 0; j--)
		{
			for (i = *pnRow - 1; i >= 0; i--)
				*(pMat1++) -= *pSubstract ;
			pSubstract++ ;
		}
	}

	void MatSubstractVectorByCol (double *pMat1, int *pnRow, int *pnCol, double *pSubstract)
	{
		int j ; 

		for (j = *pnCol - 1; j >= 0; j--)
		{
			VectorSubstractVector (pMat1, pnRow, pSubstract) ;
			pMat1 += *pnRow ;
		}
	}


	void MatSubMat (double *pMat1, double *pMat2, int *pnRow, int *pnCol)
	{
		int i; 
		for (i = *pnRow * *pnCol - 1; i >= 0; i--)
			*(pMat1++) -= *(pMat2++) ;
	}

	void VectorMult (double *x, int *pn, double dFactor)
	{
		if (dFactor == 1.0)
			return ;
		int i ;
		for (i = *pn - 1; i >=0; i--)
			x[i] *= dFactor ;
	}

	void VectorAbs (double *x, int *pn) 
	{
		int i ;
		for (i = *pn - 1; i >=0; i--)
			x[i] = ::fabs (x[i]);
	}

	void VectorSign (double *pdData, int n, double *pSgn)
	{
		int i ;
		for (i = n - 1; i >=0; i--)
		{
			if (pdData[i] > 0)
				pSgn[i] = 1 ;
			else if (pdData[i] < 0)
				pSgn[i] = -1 ;
			else
				pSgn[i] = 0 ;
		}
	}

	void InvertVector (double *pdData, int *pnLength)
	{
		int i ;
		for (i = *pnLength-1; i >= 0; i--)
		{
			*pdData = 1 / *pdData ;
			pdData++ ;
		}
	}

	void convolve (double *a, int *na, double *b, int *nb, double *ab)
	{
		int i, j, nab = *na + *nb - 1 ;
		for (i = 0; i < nab; i++)
			ab[i] = 0.0 ;
		for (i = 0; i < *na; i++)
			for (j = 0; j < *nb; j++)
				ab [i+j] += a[i] * b[j] * 100 ;
	}

	void matmult (double *pMat1, int *pnRow1, int *pnCol1, double *pMat2, int *pnCol2, double *pMatRet)
	{			
				//	verified..
		int nRow1 = *pnRow1 ;
		int nCol1 = *pnCol1 ;
		int nCol2 = *pnCol2 ;

		int i,j,h ;
		for (i = 0; i < nRow1; i++)
			for (j = 0; j  < nCol2; j++)
			{
				double *pdCurPoint = pMatRet + i + j * nRow1 ;
				*pdCurPoint = 0 ;
				for (h = 0; h < nCol1; h++)
					*pdCurPoint += pMat1 [i + h * nRow1] * pMat2 [h + j * nCol1] ;
			}
	}

	void matmult__t (double *pMat1, int *pnRow1, int *pnCol1, double *pMat2, int *pnCol2, double *pMatRet)
	{
				//	verified..			//	entpricht pMatRet = pMat1 %*% t(pMat2)
		int nRow1 = *pnRow1 ;
		int nCol1 = *pnCol1 ;
		int nCol2 = *pnCol2 ;

		int i,j,h ;
		for (i = 0; i < nRow1; i++)
			for (j = 0; j  < nCol2; j++)
			{
				double *pdCurPoint = pMatRet + i + j * nRow1 ;
				*pdCurPoint = 0 ;
				for (h = 0; h < nCol1; h++)
					*pdCurPoint += pMat1 [i + h * nRow1] * pMat2 [j + h * nCol2] ;
			}
	}

	void matmult__t_otherSize (double *pMat1, int *pnRow1, int *pnCol1, double *pMat2, int *pnCol2, double *pMatRet, int nReal2ndLength)
	{
				//	verified..			//	entpricht pMatRet = pMat1 %*% t(pMat2)
		int nRow1 = *pnRow1 ;
		int nCol1 = *pnCol1 ;
		int nCol2 = *pnCol2 ;

		int i,j,h ;
		for (i = 0; i < nRow1; i++)
			for (j = 0; j  < nCol2; j++)
			{
				double *pdCurPoint = pMatRet + i + j * nRow1 ;
				*pdCurPoint = 0 ;
				for (h = 0; h < nCol1; h++)
					*pdCurPoint += pMat1 [i + h * nRow1] * pMat2 [j + h * nReal2ndLength] ;
			}
	}

	void matmult_t_ (double *pMat1, int *pnRow1, int *pnCol1, double *pMat2, int *pnCol2, double *pMatRet)
	{
				//	verified..			//	entpricht pMatRet = t(pMat1) %*% pMat2
		int nRow1 = *pnRow1 ;
		int nCol1 = *pnCol1 ;
		int nCol2 = *pnCol2 ;

		int i,j,h ;
		for (i = 0; i < nRow1; i++)
			for (j = 0; j  < nCol2; j++)
			{
				double *pdCurPoint = pMatRet + i + j * nRow1 ;
				*pdCurPoint = 0 ;
				for (h = 0; h < nCol1; h++)
					*pdCurPoint += pMat1 [h + i * nCol1] * pMat2 [h + j * nCol1] ;
			}
	}

	void SetColOrder (double *pMat, int *pnRow, int *pnCol, int *pnOrder, double *pOrdered)
	{
		int i ;
		int nRow = *pnRow ;
		for (i = *pnCol - 1; i >= 0; i--)
			memcpy (pOrdered +  i * nRow, pMat + pnOrder[i] * nRow, sizeof (double) * nRow) ;
	}

	void SetVectorOrder (double *pdVector, int *pnSize, int *pnOrder, double *pdOrdered)
	{
		int i ;
		for (i = *pnSize - 1; i >= 0; i--)
			pdOrdered[pnOrder[i]] = pdVector[i] ;
	}

	void sd (double *pData, int *pn, double *pdSd)
	{
		int i ;
		double dMean ;

		mean (pData, pn, &dMean) ;

		double dCurDistance ;
		double dSum = 0 ;
		for (i = *pn - 1; i >= 0; i--)
		{
			dCurDistance = *pData++ - dMean ;
			dSum += dCurDistance * dCurDistance ;

		}
		*pdSd = ::sqrt (dSum / (*pn - 1)) ;
	}

	void mean (double *pData, int *npLength, double *pMean)
	{
		double dSum = 0 ;
		int i ;
		for (i = *npLength - 1; i >= 0; i--)
			dSum += pData [i] ;
		*pMean = dSum / *npLength ;
	}

	void ColMeans (double *pData, int *pnRows, int *pnCols, double *pMean)
	{
		int i ;
		for (i = *pnCols - 1; i >= 0; i--)
			mean (pData + i * *pnRows, pnRows, pMean + i) ;
	}

	void colsd (double *pdMat, int *npRow, int *npCol, double *pdSd)
	{
				//	verified..
		int nRow = *npRow ;
		int j ;
		for (j = *npCol - 1; j >= 0; j--)
			sd (pdMat + (j * nRow), npRow, pdSd + j) ;
	}

	void fillDiag (double *pdMat, int *npRow, int *npCol)
	{
		int i,j ;
		int nRow = *npRow;
		for (i = *npCol - 1; i >= 0; i--)
			for (j = nRow - 1; j >= 0; j--)
				if (i != j)
					pdMat [i * nRow + j] = 0 ;
				else
					pdMat [i * nRow + j] = 1 ;
	}

	void fillDiag_Value (double *pdMat, int *npRow, int *npCol, double *pdValue)
	{
		int i,j ;
		int nRow = *npRow;
		for (i = *npCol - 1; i >= 0; i--)
			for (j = nRow - 1; j >= 0; j--)
				if (i != j)
					pdMat [i * nRow + j] = 0 ;
				else
					pdMat [i * nRow + j] = *pdValue ;
	}

	void fill (double *pdData, int nLength, double dValue)
	{
		int i ;
		for (i = nLength - 1; i >= 0; i--)
			pdData[i] = dValue ;
	}

	void median (double *pData, int *npLength, double *pMedian, int bData_ReadOnly)
	{			//	verified...
		int n = *npLength ;
		if (!n) 
		{
			*pMedian = R_NaN ;
			return ;
		}
		int nHalf = (n + 1)/2 ;

		double *pData_Copy ;

		if (bData_ReadOnly)
		{
			pData_Copy = (double *) malloc (sizeof (double) * n) ;
			memcpy (pData_Copy, pData, sizeof (double) * n) ;
		}
		else
			pData_Copy = pData ;

		if (n % 2)
		{
			rPsort(pData_Copy, n, nHalf-1);
			*pMedian = pData_Copy[nHalf - 1] ;
		}
		else
		{
			double dTemp ;
			rPsort(pData_Copy, n, nHalf - 1);
			dTemp = pData_Copy[nHalf - 1] ;
			rPsort(pData_Copy, n, nHalf);

			*pMedian = (dTemp + pData_Copy[nHalf]) / 2 ;
		}


/*		R_rsort (pData_Copy, n) ;

		if (n % 2)
			*pMedian = pData_Copy[nHalf - 1] ;
		else
			*pMedian = (pData_Copy[nHalf - 1] + pData_Copy[nHalf]) / 2 ;
*/
		if (bData_ReadOnly)
			free (pData_Copy) ;
	}

	double max (double *pdData, int nSize)
	{
		double dRet = *pdData ;
		int i ;
		for (i = nSize - 1; i; i--)
			if (dRet < pdData[i])
				dRet = pdData[i] ;
		return dRet ;
	}

	double median (double *pdData, int nSize)
	{
		double *pdWork = (double *) malloc (nSize * sizeof(double)) ;
		memcpy (pdWork, pdData, nSize * sizeof (double)) ;
		double dMedian = median_raw (pdWork, nSize) ;
		free (pdWork) ;
		return dMedian ;
	}

	double median (double *pdData, double *pdWork, int nSize)
	{
		memcpy (pdWork, pdData, nSize * sizeof (double)) ;
		return median_raw (pdWork, nSize) ;
	}

	double median_raw (double *pdWork, int nSize)
	{
		int nHalf = nSize >> 1 ;
		rPsort(pdWork, nSize, nHalf) ;
		if (nSize % 2)
			return pdWork [nHalf] ;
		//rPsort(pdWork, nSize, nHalf + 1);
		return (pdWork [nHalf] + max (pdWork, nHalf)) / 2 ;
		//return (pdWork [nHalf] + pdWork [nHalf + 1]) / 2 ;
	}

	void ColMedian (double *pData, int *pnRows, int *pnCols, double *pMedians, int bData_ReadOnly)
	{
		int i ;
		int nRows = *pnRows ;
		for (i = *pnCols - 1; i >= 0; i--)
			median (pData + i * nRows, pnRows, pMedians + i, bData_ReadOnly) ;
	}

	struct OrderStruct
	{
		double	m_dValue ;
		int		m_nOrder ;
	} ;

	int compare (struct OrderStruct *elem1, struct OrderStruct *elem2 )
	{
		if (elem1->m_dValue < elem2->m_dValue)
			return -1 ;
		if (elem1->m_dValue > elem2->m_dValue)
			return 1 ;
		return 0 ;
	}

	void order  (double *pData, int *npLength, int*pnOrder)
	{	//	is actually done by R as well... ;)
		struct OrderStruct *pOrder = (struct OrderStruct *) malloc (*npLength * sizeof (struct OrderStruct)) ;
		int i ;
		for (i = *npLength - 1; i >= 0; i--)
		{
			pOrder[i].m_dValue = pData[i] ;
			pOrder[i].m_nOrder = i ;
		}

		qsort (pOrder, *npLength, sizeof (struct OrderStruct), (int ( *)(const void *,const void *)) compare) ;

		for (i = *npLength - 1; i >= 0; i--)
			pnOrder[i] = pOrder[i].m_nOrder ;

		free (pOrder) ;
	}

	void order_decr  (double *pData, int *npLength, int*pnOrder)
	{
		struct OrderStruct *pOrder = (struct OrderStruct *) malloc (*npLength * sizeof (struct OrderStruct)) ;
		int i ;
		for (i = *npLength - 1; i >= 0; i--)
		{
			pOrder[i].m_dValue = pData[i] ;
			pOrder[i].m_nOrder = i ;
		}

		qsort (pOrder, *npLength, sizeof (struct OrderStruct), (int ( *)(const void *,const void *)) compare) ;

		for (i = *npLength - 1; i >= 0; i--)
			*pnOrder++ = pOrder[i].m_nOrder ;

		free (pOrder) ;
	}

	void mad (double *x, int *npLength, double *pdMad) 
	{			// verified...
//		std - mad: als center wird der median angenommen, die Constante ist fix 1.4826, nas werden nicht gelöscht. high & low sind auf false gesetzt1
		int n = *npLength ;

		double dCenter ;

		double *px_Copy = (double *) malloc (sizeof (double) * n) ;
		memcpy (px_Copy, x, sizeof (double) * n) ;

		median (px_Copy, npLength, &dCenter, 0) ;

		VectorAdd (px_Copy, npLength, -dCenter) ;
		VectorAbs (px_Copy, npLength) ;

		

		median (px_Copy, npLength, pdMad, 0) ;

		free (px_Copy) ;

		*pdMad *= 1.4826 ;
	}

	void rowmad (double *pdMat, int *pnRow, int *pnCol, double *pdMad)
	{
		double *pdtMat = new double [*pnRow * *pnCol] ;
		t (pdMat, pnRow, pnCol, pdtMat) ;
		colmad (pdtMat, pnCol, pnRow, pdMad) ;
		delete [] pdtMat ;
	}

	void rowsd (double *pdMat, int *pnRow, int *pnCol, double *pdSd)
	{
		double *pdtMat = new double [*pnRow * *pnCol] ;
		t (pdMat, pnRow, pnCol, pdtMat) ;
		colsd (pdtMat, pnCol, pnRow, pdSd) ;
		delete [] pdtMat ;
	}

	void rowdev (double *pdMat, int *pnRow, int *pnCol, double *pdDev, dev_Type *pfDev)
	{
		double *pdtMat = new double [*pnRow * *pnCol] ;
		t (pdMat, pnRow, pnCol, pdtMat) ;
		coldev (pdtMat, pnCol, pnRow, pdDev, pfDev) ;
		delete [] pdtMat ;
	}

	void coldev (double *pdMat, int *npRow, int *npCol, double *pdDev, dev_Type *pfDev)
	{
				//	verified..
		int nRow = *npRow ;
		int j ;
		for (j = *npCol - 1; j >= 0; j--)
			pfDev (pdMat + (j * nRow), npRow, pdDev + j) ;
	}

	void colmad (double *pdMat, int *npRow, int *npCol, double *pdSd)
	{
				//	verified..
		int nRow = *npRow ;
		int j ;
		for (j = *npCol - 1; j >= 0; j--)
			mad (pdMat + (j * nRow), npRow, pdSd + j) ;
	}

	void t (double *pSrc, int *pnRow, int *pnCol, double *pDest)
	{			// verified...
		int i, j ;
		int nRow = *pnRow, nCol = *pnCol ;

		for (i = nRow - 1; i >= 0; i--)
			for (j = nCol - 1; j >= 0; j--)
				pDest [j + i * nCol] = pSrc [i + j * nRow] ;
	}

	void wheremax (double *pData, int *pnSize, int *pnIdx)
	{
		*pnIdx = 0 ;
		double dMax = pData[0] ;
		int i ;
		for (i = *pnSize - 1; i > 0; i--)
			if (pData[i] > dMax)
			{
				dMax = pData[i] ;
				*pnIdx = i ;
			}
	}


	void cbind (double *pMat1, int *pnRow, int *pnCol1, double *pMat2, int *pnCol2, double *pDestMat)
	{
		int i ;
		int nRow = *pnRow ;
		for (i = *pnCol1 - 1; i >= 0; i--)
			memcpy (pDestMat + i * nRow, pMat1 + i * nRow, sizeof(double) * nRow) ;

		int nOffset = *pnCol1 ;
		for (i = *pnCol2 - 1; i >= 0; i--)
			memcpy (pDestMat + (i + nOffset) * nRow, pMat2 + i * nRow, sizeof(double) * nRow) ;
	}

	void matcpy (double *pMatDest, double *pMatSource, int nSize)
	{
		memcpy (pMatDest, pMatSource, nSize * sizeof(double)) ;
	}

	void center (double *pdData, int *pnRow, int *pnCol)
	{
		int i ;
		double dMean ;
		int nRow = *pnRow ;
		for (i = *pnCol - 1; i >= 0; i--)
		{
			mean (pdData + i * nRow, pnRow, &dMean) ;
			VectorAdd (pdData + i * nRow, pnRow, dMean) ;
		}
	}

	void ColSums (double *pdData, int *pnRow, int *pnCol, double *pColSums)
	{
		int i,j ;
		for (i = *pnCol - 1; i >= 0; i--)
		{
			*pColSums = 0 ;
			for (j = *pnRow - 1; j >= 0; j--)
				*pColSums += *pdData++ ;
			pColSums++ ;
		}
	}

	void GetRow (double *pdData, int *pnRow, int *pnCol, int *pnRowIdx, double *pDest)
	{
		int nRow = *pnRow ;
		int i ; 
		pdData += *pnRowIdx ;
		for (i = *pnCol - 1; i >= 0; i--)
		{
			*pDest++ = *pdData ;
			pdData+= nRow ;
		}

	}

	double VectorMin (double *pdData, int *pnLength)
	{
		double dRet = *pdData ;
		int i ;
		for (i = *pnLength - 1; i; i--)
			if (pnLength [i] < dRet)
				dRet = pnLength [i] ;
		return dRet ;

	}

	int	which_max (double *pdData, int *pnLength)
	{
		int nMax = 0 ;
		double dMax = pdData [0] ;
		int i ;
		for (i = *pnLength - 1; i > 0; i--)
			if (pdData[i] > dMax)
			{
				dMax = pdData[i] ;
				nMax  = i ;
			}
		return nMax ;
	}

	void myrand (int *pn, double *pdR0, double *pdRet)
	{
		double dRn = *pdR0 ;
		int i ;
		double foo ;
//		TRACE ("%0.20f", dRn) ;
		for (i = *pn - 1; i >= 0; i--)
			*pdRet++ = dRn = modf (9821 * dRn + (double) 0.21327, &foo) ;
	}

	void VectorSqrt (double *pdData, int *pnLength, double *pdSqrt)
	{
		int i ;
		for (i = *pnLength - 1; i >= 0; i--)
			pdSqrt [i] = ::sqrt (pdData[i]) ;
	}

	void RepData (double *pdData, int nLength, int nCount, double *pdDest)
	{
		int i ;
		for (i = 0; i < nCount; i++)
		{
			memcpy (pdDest + i * nLength, pdData, nLength * sizeof (double)) ;
		}
	}

	int GetRows (double *pdX, int n, int p, int *pnIdx, double *pdOut)
	{
		int nOutN = 0 ;
		int i, j ;
		for (i = 0; i < n; i++)
			if (pnIdx[i])
				nOutN ++ ;

		for (i = 0; i < n; i++)
			if (pnIdx [i])
			{
				for (j = 0; j < p; j++)
					pdOut [j * nOutN] = pdX [i + j * n] ;
				pdOut++ ;
			}
		return nOutN ;
	}

	void GetVectorIdx (double *pdX, int n, int *pnIdx, double *pdOut)
	{
		int i ;
		for (i = n; i ; i--)
		{
			if (*pnIdx++)
				*pdOut++ = *pdX ;
			pdX++ ;
		}
	}

	dev_Type *GetDevFunction (int nDev)
	{
		switch (nDev)
		{
		case 0: return mad ;
		case 1: return sd ;
		case 2:	return Qn ;
		}
		return NULL ;
	}


	void MatSetMat (double *pdBig, double *pdSmall, int nNB, /*int nPB,*/ int nNS, int nPS, int nOffsetN, int nOffsetP)
	{				//	kopiert eine kleine nNS x nPS matrix in eine größere nNB x nPB matrix, wobei der Paramater nPB eigentlich nur für checks verwendet werden kann -deshalb wird er hier weggelassen 
		int i ;
		pdBig += nNB * nOffsetP + nOffsetN ;

		for (i = nPS - 1; i >= 0; i--)
		{
			memcpy (pdBig, pdSmall, nNS * sizeof (double)) ;
			pdBig += nNB ;
			pdSmall += nNS ;
		}
	}

	void SwapPtrs (double **pptr1, double **pptr2)
	{
		double *ptr3 ;
		ptr3 = *pptr1 ;
		*pptr1 = *pptr2 ;
		*pptr2 = ptr3 ;
	}
