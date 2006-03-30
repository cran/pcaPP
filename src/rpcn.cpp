#include "math.h"
#include "R.h"
#include "rsubst.h"
#include "fastpca.h"

void rpcn (double *pdData, int *pnParams, double *pdScores, double *pdVeig, double *pdLambda)
{								//	validated..

	int &n = pnParams[0] ;
	int &p = pnParams[1] ;
	int &k = pnParams[2] ;
	int &nScal = pnParams[3] ;
	int &nScores = pnParams[4] ;
	int &nRealN = pnParams[5] ;


	dev_Type *pfDev = GetDevFunction(nScal) ;


	double *pdy = pdData ;

	double *pdA = new double [n * p] ;		//	n x p

	double *pdY = new double [n * n] ;		//	n x n

	double *pdHlp = new double [n] ;		//	n
	int	   *pnHlpTF = new int [n] ;			//	n

	double *pdPcol = pdHlp ; //new double [n] ;		//	n

	double *pdScoreVec = pdHlp ;//new double [n] ;	//

	double *pdScore_x_Veig = pdA ; //new double [n * p] ;	//	n x p

	int nOne = 1 ;

	int i ;
	int nShortN ;

	double *pdCurLambda ;
	double *pdCurVeig ;

	for (i = 0; i < k; i++)
	{
		pdCurLambda = pdLambda  + i ;
		pdCurVeig = pdVeig + i * p ;

		//	hlp = colSums (t(y)^2)										//	in: pdHlp
		calcrowsumsq (pdy, n, p, pdHlp) ;
		VectorLargerThan (pdHlp, n, pnHlpTF, 0.0000000000000001) ;
		calcvectorsqrt (pdHlp, n, pdHlp) ;

		//	A = y[hlpTF,] / sqrt (hlp[lpTF])
		GetVectorIdx (pdHlp, n, pnHlpTF, pdHlp) ;

		nShortN = GetRows (pdy, n, p, pnHlpTF, pdA) ;														//	in: pdA
		MatDivVectorByCol (pdA, &nShortN, &p, pdHlp) ;					//	out: pdHlp

		matmult__t_otherSize (pdA, &nShortN, &p, pdy, &nRealN, pdY, n) ;// (nShortN x p) * (p * nRealN)

			//	pcol=feval(s,Y');
/*		if (nScal == 0)
			rowsd (pdY, &nShortN, &nRealN, pdPcol) ;					//	in: pdPcol
		else if (nScal == 1)
			rowmad (pdY, &nShortN, &nRealN, pdPcol) ;*/
		rowdev (pdY, &nShortN, &nRealN, pdPcol, pfDev) ;

			//	[lambdastar,istar]=max(pcol);
		int nIStar = which_max (pdPcol, &nShortN) ;
		*pdCurLambda = pdPcol[nIStar] ;									//	out: pdPcol

		GetRow (pdA, &nShortN, &p, &nIStar, pdCurVeig) ;													//	out: pdA

																		//	in: pdScoreVec
		matmult (pdy, &n, &p, pdCurVeig, &nOne, pdScoreVec) ;

		if (nScores)
			matcpy (pdScores + nRealN * i, pdScoreVec, nRealN) ;

		if (i < k - 1)																					//	in: pdScore_x_Veig
		{
			matmult (pdScoreVec, &n, &nOne, pdCurVeig, &p, pdScore_x_Veig) ;
			MatSubMat (pdy, pdScore_x_Veig, &n, &p) ;
		}																//	out: pdScoreVec				//	out: pdScore_x_Veig
	}

//	VectorMultVector (pdLambda, &k, pdLambda) ;

	delete [] pdScore_x_Veig ;

	delete [] pdHlp ;
	delete [] pnHlpTF ;

	delete [] pdY ;
}
