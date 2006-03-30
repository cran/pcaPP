#include "math.h"
#include "R.h"
#include "rsubst.h"
#include "fastpca.h"

void rpcnup (double *pdData, int *pnParams, double *pdScores, double *pdVeig, double *pdLambda)
{
	int &n = pnParams[0] ;
	int &p = pnParams[1] ;
	int &k = pnParams[2] ;
	int &nScal = pnParams[3] ;
	int &nScores = pnParams[4] ;

	int nMaxit = pnParams[5] ;
	int nMaxhalf = pnParams[6] ;

	int nRealN = pnParams[7] ;

	dev_Type *pfDev = GetDevFunction(nScal) ;

	double *pdy = pdData ;//new double [n * p] ;		//	n x p

	double *pdA = new double [n * p] ;		//	n x p

	double *pdY = new double [n * n] ;		//	n x n
	double *pdPcol = new double [n] ;		//	n

	double *pdVH = new double [p] ;			//	p

	double *pdTScores = new double [n] ;	//	n

	double *pdHlp = new double [n] ;		//	n
	int	   *pnHlpTF = new int [n] ;			//	n

	double *pdScore_x_Veig = new double [n * p] ;	//	n x p		//	use pdA ?!

	double *pdCurRowA = new double [p] ;

	int nOne = 1 ;

	int i, m, kk ;

	int nShortN ;

	double dNewObj ;

	double *pdCurLambda ;
	double *pdCurVeig ;
//	double *pdCurScores ;
//	if (!nScores)
//		pdCurScores = new double [n] ;

	double *pdScoreVec = new double [n] ;

	double *pdScoreSign = new double [n] ;	//	n

	for (i = 0; i < k; i++)
	{
		pdCurLambda = pdLambda  + i ;
		pdCurVeig = pdVeig + i * p ;

		//	hlp = colSums (t(y)^2)
		calcrowsumsq (pdy, n, p, pdHlp) ;
		VectorLargerThan (pdHlp, n, pnHlpTF, 0.0000000000000001) ;

//pnHlpTF[0] = 0 ;
		calcvectorsqrt (pdHlp, n, pdHlp) ;

		//	A = y[hlpTF,] / sqrt (hlp[hlpTF])
		GetVectorIdx (pdHlp, n, pnHlpTF, pdHlp) ;

		nShortN = GetRows (pdy, n, p, pnHlpTF, pdA) ;														//	in: pdA
		MatDivVectorByCol (pdA, &nShortN, &p, pdHlp) ;			//	out: pdHlp

		matmult__t_otherSize (pdA, &nShortN, &p, pdy, &nRealN, pdY, n) ;// (nShortN x p) * (p * nRealN)

//		matmult__t (pdA, &nShortN, &p, pdy, &n, pdY) ;// (nShortN x p) * (p * n)

		rowdev (pdY, &nShortN, &nRealN, pdPcol, pfDev) ;
//		rowdev (pdY, &nShortN, &n, pdPcol, pfDev) ;

			//	[lambdastar,istar]=max(pcol);
		int nIStar = which_max (pdPcol, &nShortN) ;
		*pdCurLambda = pdPcol[nIStar] ;

		GetRow (pdA, &nShortN, &p, &nIStar, pdCurVeig) ;

//		if (nScores)
//			pdCurScores = pdScores + i * n ;
//		//	pdCurScores = pdScores + nRealN * i

		matmult (pdy, &n, &p, pdCurVeig, &nOne, pdScoreVec) ;

		///	updating procedure

		for (m = 0; m < nMaxit; m++)
		{
			VectorSign (pdScoreVec, n, pdScoreSign) ;

//			t(pdA, &n, &p, pdAt) ;
//			matmult (pdAt, &p, &n, pdScoreSign, &nOne, pdVH) ;

			GetVectorIdx (pdScoreSign, n, pnHlpTF, pdScoreSign) ;

			matmult_t_ (pdA, &p, &nShortN, pdScoreSign, &nOne, pdVH) ;		//	p x nShortN * nShortN x 1

			//vh=vh/norm(vh);
			//vh=vh/sqrt(sum(vh^2))
			VectorMult (pdVH, &p, 1 / (double)nShortN) ;

			VectorMult (pdVH, &p, 1 / calcsqrtsumsqr (pdVH, p)) ;

			//	vh = (t(A)%*%sign(scorevec[hlpTF]))/n
			//	vh = vh / sqrt(sum(vh^2))

			//tscores=y*vh;
			matmult (pdy, &n, &p, pdVH, &nOne, pdTScores) ;

			//newobj=feval(s,tscores);

			pfDev (pdTScores, &n, &dNewObj) ;

			kk=0;

			while (dNewObj < *pdCurLambda  && kk < nMaxhalf)
			{
//				vh=(vh+vhelp)/2;
				VectorAddVector (pdVH, &p, pdCurVeig) ;
				VectorMult (pdVH, &p, 0.5) ;
//				vh=vh/norm(vh);
				VectorMult (pdVH, &p, 1 / calcsqrtsumsqr (pdVH, p)) ;

				//tscores=y*vh;
				matmult (pdy, &n, &p, pdVH, &nOne, pdTScores) ;

				//newobj = feval (s, tscores);

				pfDev (pdTScores, &n, &dNewObj) ;

				kk++ ;
			}

			if (dNewObj > *pdCurLambda)
			{
				matcpy (pdScoreVec, pdTScores, n) ;
			//	matcpy (pdCurScores, pdTScores, nRealN) ;
				matcpy (pdCurVeig, pdVH, p) ;
				*pdCurLambda = dNewObj;
			}
			else
				break ;
		}

		if (nScores)
			matcpy (pdScores + i * nRealN, pdScoreVec, nRealN) ;

		if (i < k - 1)
		{
			
//		scorevec %*% A[istar,]
	//		GetRow (pdA, &nShortN, &p, &nIStar, pdCurRowA) ;
	//		matmult (pdScoreVec, &n, &nOne, pdCurRowA, &p, pdScore_x_Veig) ;
//		scorevec %*% t(vhelp)	//	da vhelp nie unabhängig von pdCurVeig verwendet wird --> selber membereich
			matmult (pdScoreVec, &n, &nOne, pdCurVeig, &p, pdScore_x_Veig) ;
//	weiter..

			MatSubMat (pdy, pdScore_x_Veig, &n, &p) ;
		}
	}

//	VectorMultVector (pdLambda, &k, pdLambda) ;

//	if (!nScores)
//		delete [] pdCurScores ;

	delete [] pdA ;
	
	delete [] pdY ;
	delete [] pdPcol ;

	delete [] pdVH ;

	delete [] pdTScores ;
	
	delete [] pdHlp ;
	delete [] pnHlpTF ;

	delete [] pdScore_x_Veig ;

	delete [] pdCurRowA ;

	delete [] pdScoreSign ;

	delete [] pdScoreVec ;
}


/*

  void rpcnup (double *pdData, int *pnParams, double *pdScores, double *pdVeig, double *pdLambda)
{
	int &n = pnParams[0] ;
	int &p = pnParams[1] ;
	int &k = pnParams[2] ;
	int &nScal = pnParams[3] ;
	int &nScores = pnParams[4] ;

	int nMaxit = pnParams[5] ;
	int nMaxhalf = pnParams[6] ;

	int nRealN = pnParams[7] ;

	dev_Type *pfDev = GetDevFunction(nScal) ;

	double *pdy = pdData ;//new double [n * p] ;		//	n x p

	double *pdA = new double [n * p] ;		//	n x p

	double *pdY = new double [n * n] ;		//	n x n
	double *pdPcol = new double [n] ;		//	n

	double *pdVH = new double [p] ;			//	p

	double *pdTScores = new double [n] ;	//	n

	double *pdHlp = new double [n] ;		//	n
	int	   *pnHlpTF = new int [n] ;			//	n

	double *pdScore_x_Veig = new double [n * p] ;	//	n x p		//	use pdA ?!

	double *pdCurRowA = new double [p] ;

	int nOne = 1 ;

	int i, m, kk ;

	int nShortN ;

	double dNewObj ;

	double *pdCurLambda ;
	double *pdCurVeig ;
	double *pdCurScores ;
	if (!nScores)
		pdCurScores = new double [n] ;

	double *pdScoreSign = new double [n] ;	//	n

	for (i = 0; i < k; i++)
	{
		pdCurLambda = pdLambda  + i ;
		pdCurVeig = pdVeig + i * p ;

		//	hlp = colSums (t(y)^2)
		calcrowsumsq (pdy, n, p, pdHlp) ;
		VectorLargerThan (pdHlp, n, pnHlpTF, 0.0000000000000001) ;

//pnHlpTF[0] = 0 ;
		calcvectorsqrt (pdHlp, n, pdHlp) ;

		//	A = y[hlpTF,] / sqrt (hlp[hlpTF])
		GetVectorIdx (pdHlp, n, pnHlpTF, pdHlp) ;

		nShortN = GetRows (pdy, n, p, pnHlpTF, pdA) ;														//	in: pdA
		MatDivVectorByCol (pdA, &nShortN, &p, pdHlp) ;			//	out: pdHlp

//		matmult__t_otherSize (pdA, &nShortN, &p, pdy, &nRealN, pdY, n) ;// (nShortN x p) * (p * nRealN)

		matmult__t (pdA, &nShortN, &p, pdy, &n, pdY) ;// (nShortN x p) * (p * n)

//		rowdev (pdY, &nShortN, &nRealN, pdPcol, pfDev) ;
		rowdev (pdY, &nShortN, &n, pdPcol, pfDev) ;

			//	[lambdastar,istar]=max(pcol);
		int nIStar = which_max (pdPcol, &nShortN) ;
		*pdCurLambda = pdPcol[nIStar] ;

		GetRow (pdA, &nShortN, &p, &nIStar, pdCurVeig) ;

		if (nScores)
			pdCurScores = pdScores + i * n ;
		//	pdCurScores = pdScores + nRealN * i

		matmult (pdy, &n, &p, pdCurVeig, &nOne, pdCurScores) ;

		///	updating procedure

		for (m = 0; m < nMaxit; m++)
		{
			VectorSign (pdCurScores, n, pdScoreSign) ;

//			t(pdA, &n, &p, pdAt) ;
//			matmult (pdAt, &p, &n, pdScoreSign, &nOne, pdVH) ;

			GetVectorIdx (pdScoreSign, n, pnHlpTF, pdScoreSign) ;

			matmult_t_ (pdA, &p, &nShortN, pdScoreSign, &nOne, pdVH) ;		//	p x nShortN * nShortN x 1

			//vh=vh/norm(vh);
			//vh=vh/sqrt(sum(vh^2))
			VectorMult (pdVH, &p, 1 / (double)nShortN) ;

			VectorMult (pdVH, &p, 1 / calcsqrtsumsqr (pdVH, p)) ;

			//	vh = (t(A)%*%sign(scorevec[hlpTF]))/n
			//	vh = vh / sqrt(sum(vh^2))

			//tscores=y*vh;
			matmult (pdy, &n, &p, pdVH, &nOne, pdTScores) ;

			//newobj=feval(s,tscores);

			pfDev (pdTScores, &n, &dNewObj) ;

			kk=0;

			while (dNewObj < *pdCurLambda  && kk < nMaxhalf)
			{
//				vh=(vh+vhelp)/2;
				VectorAddVector (pdVH, &p, pdCurVeig) ;
				VectorMult (pdVH, &p, 0.5) ;
//				vh=vh/norm(vh);
				VectorMult (pdVH, &p, 1 / calcsqrtsumsqr (pdVH, p)) ;

				//tscores=y*vh;
				matmult (pdy, &n, &p, pdVH, &nOne, pdTScores) ;

				//newobj=feval(s,tscores);

//				if (nScal == 0)
//					sd (pdTScores, &n, &dNewObj) ;
//				else if (nScal == 1)
//					mad (pdTScores, &n, &dNewObj) ;
//
				pfDev (pdTScores, &n, &dNewObj) ;

				kk++ ;
			}

			if (dNewObj > *pdCurLambda)
			{
				matcpy (pdCurScores, pdTScores, n) ;
				matcpy (pdCurVeig, pdVH, p) ;
				*pdCurLambda = dNewObj;
			}
			else
				break ;
		}

		if (i < k - 1)
		{
			GetRow (pdA, &nShortN, &p, &nIStar, pdCurRowA) ;
			matmult (pdCurScores, &n, &nOne, pdCurRowA, &p, pdScore_x_Veig) ;
			MatSubMat (pdy, pdScore_x_Veig, &n, &p) ;
		}
	}

	VectorMultVector (pdLambda, &k, pdLambda) ;

	if (!nScores)
		delete [] pdCurScores ;

	delete [] pdA ;
	
	delete [] pdY ;
	delete [] pdPcol ;

	delete [] pdVH ;

	delete [] pdTScores ;
	
	delete [] pdHlp ;
	delete [] pnHlpTF ;

	delete [] pdScore_x_Veig ;

	delete [] pdCurRowA ;

	delete [] pdScoreSign ;
}


*/
