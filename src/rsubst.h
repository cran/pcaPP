#ifdef MSDEBUGTEST
	#pragma comment(lib, "Rdll.lib") 
#endif // #ifdef MSDEBUGTEST


#define NEWDARRAY(COUNT) ((double *) malloc ((COUNT) * sizeof (double)))
#define NEWNARRAY(COUNT) ((int *) malloc ((COUNT) * sizeof (int)))

#define min(a,b)            (((a) < (b)) ? (a) : (b))

typedef int                 BOOL;

#define FALSE 0
#define TRUE 1
#ifndef NULL
	#define NULL 0
#endif


	typedef void (dev_Type)(double *pData, int *pn, double *pdRet) ;

	//	R-Functionst
//	void VectorAbs (double *x, int *n)  ;

	// _NOT_ R-Functions
	void median (double *pData, int *npLength, double *pMedian, int bData_ReadOnly) ;

	void VectorAdd (double *x, int *pnSize, double dAdd) ;

	void VectorAddVector (double *x, int *pn, double *pdAdd) ;
	void VectorSubstractVector (double *x, int *pn, double *pdSubstract) ;
	void VectorMultVector (double *x, int *pn, double *pdMult) ;
	void VectorDivVector (double *x, int *pn, double *pdDiv) ;
	void VectorDivVectorIdx (double *x, int *pn, double *pdDiv, int *pnIdx) ;
	void MatDivVectorByColIdx (double *pdData, int n, int p, double *pdDiv, int *pdIdx) ;
	void MatAddMat (double *pMat1, double *pMat2, int *pnRow, int *pnCol) ;
	void MatSubMat (double *pMat1, double *pMat2, int *pnRow, int *pnCol) ;


	void MatAddVectorByRow (double *pMat1, int *pnRow, int *pnCol, double *pAdd) ;
	void MatMultVectorByCol (double *pdData, int *pnRow, int *pnCol, double *pMult) ;
	void MatSubstractVectorByCol (double *pMat1, int *pnRow, int *pnCol, double *pSubstract) ;
	void MatDivVectorByCol (double *pdData, int *pnRow, int *pnCol, double *pDiv) ;
	void MatSubstractVectorByRow (double *pMat1, int *pnRow, int *pnCol, double *pSubstract) ;

	void VectorMult (double *x, int *pn, double dFactor) ;
	void InvertVector (double *pdData, int *pnLength) ;
	void VectorAbs (double *x, int *pn) ;
	void VectorSign (double *pdData, int n, double *pSgn) ;

	void convolve (double *a, int *na, double *b, int *nb, double *ab) ;
	void matmult (double *pMat1, int *pnRow1, int *pnCol1, double *pMat2, int *pnCol2, double *pMatRet) ;
	void matmult__t (double *pMat1, int *pnRow1, int *pnCol1, double *pMat2, int *pnCol2, double *pMatRet) ;
	void matmult_t_ (double *pMat1, int *pnRow1, int *pnCol1, double *pMat2, int *pnCol2, double *pMatRet) ;
	void matmult__t_otherSize (double *pMat1, int *pnRow1, int *pnCol1, double *pMat2, int *pnCol2, double *pMatRet, int nReal2ndLength) ;


	void mean (double *pData, int *npLength, double *pMean) ;
	void ColMeans (double *pData, int *pnRows, int *pnCols, double *pMean) ;

	void SetColOrder (double *pMat, int *pnRow, int *pnCol, int *pnOrder, double *pOrdered) ;
	void SetVectorOrder (double *pdVector, int *pnSize, int *pnOrder, double *pdOrdered) ;
	void sd (double *pData, int *pn, double *pdSd) ;
	void colsd (double *pdMat, int *npRow, int *npCol, double *pdSd) ;
	void rowsd (double *pdMat, int *pnRow, int *pnCol, double *pdSd) ;
	void fillDiag (double *pdMat, int *npRow, int *npCol) ;
	void fillDiag_Value (double *pdMat, int *npRow, int *npCol, double *pdValue) ;
	void fill (double *pdData, int nLength, double dValue) ;
	void median (double *pData, int *npLength, double *pMedian, int bData_ReadOnly) ;
	void ColMedian (double *pData, int *npRows, int *npCols, double *pMedians, int bData_ReadOnly) ;

	double median (double *pdData, int nSize) ;
	double median (double *pdData, double *pdWork, int nSize) ;
	double median_raw (double *pdWork, int nSize) ;

	int compare (struct OrderStruct *elem1, struct OrderStruct *elem2 ) ;

	void order  (double *pData, int *npLength, int*pnOrder) ;
	void order_decr  (double *pData, int *npLength, int*pnOrder) ;
	void mad (double *x, int *npLength, double *pdMad) ;
	void colmad (double *pdMat, int *npRow, int *npCol, double *pdSd) ;

	void coldev (double *pdMat, int *npRow, int *npCol, double *pdSd, dev_Type *pfDev) ;
	void rowdev (double *pdMat, int *pnRow, int *pnCol, double *pdDev, dev_Type *pfDev) ;

	void rowmad (double *pdMat, int *pnRow, int *pnCol, double *pdMad) ;
	void t (double *pSrc, int *pnRow, int *pnCol, double *pDest) ;
	void wheremax (double *pData, int *pnSize, int *pnIdx) ;
	void cbind (double *pMat1, int *pnRow, int *pnCol1, double *pMat2, int *pnCol2, double *pDestMat) ;
	void matcpy (double *pMatDest, double *pMatSource, int nSize) ;
	void center (double *pdData, int *pnRow, int *pnCol) ;

	void ColSums (double *pdData, int *pnRow, int *pnCol, double *pColSums) ;

	void GetRow (double *pdData, int *pnRow, int *pnCol, int *pnRowIdx, double *pDest) ;

	double VectorMin (double *pdData, int *pnLength) ;

	int	which_max (double *pdData, int *pnLength) ;
	void RepData (double *pdData, int nLength, int nCount, double *pdDest) ;

	void VectorSqrt (double *pdData, int *pnLength, double *pdSqrt) ;


	void calcvectorsqrt (double *pdData, int n, double *pdSqrt) ;
	double calcsqrtsumsqr (double *pData, int nCount) ;
	double calcsumsqrVectorDiff (double *pD1, double *pD2, int nCount) ;
	double calcrowsumsq (double *pdData, int n, int p, double *pdRowSumSq) ;
	int VectorLargerThan (double *pdData, int n, int *pnTF, double pMargin) ;
	int GetRows (double *pdX, int n, int p, int *pnIdx, double *pdOut) ;
	void GetVectorIdx (double *pdX, int n, int *pnIdx, double *pdOut) ;

	void sd (double *pData, int *pn, double *pdSd) ;


	dev_Type *GetDevFunction (int nDev) ;


	void gridplane (double *pdY, int *npYSize, int *pm, dev_Type *pfDev, double *pdObjMax, double *pdAlphaMax) ;
	void gridplane2 (double *pdY, int *npYSize, double *pa, int *pm, double *pdDiv, dev_Type *pfDev, double *pdObjMax, double *pdAlphaMax) ;

	void MatSetMat (double *pdBig, double *pdSmall, int nNB, /*int nPB,*/ int nNS, int nPS, int nOffsetN, int nOffsetP) ;
	void SwapPtrs (double **pptr1, double **pptr2) ;
