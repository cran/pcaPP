#include "math.h"
#include "R.h"
#include "rsubst.h"
#include "fastpca.h"

	void InitD (double *pD, int n, double dVal)
	{
		int i ;
		for (i = 0; i < n; i++)
			pD[i] = dVal ;
	}

	void InitN (int *pN, int n, int nVal)
	{
		int i ;
		for (i = 0; i < n; i++)
			pN[i] = nVal ;
	}



/* Time-efficient algorithm for the scale estimator :

  Qn = dn * 2.2219 * (|x_i - x_j|; i<j)_(k)

  Parameters of the function Qn :
     x : real array containing the observations
     n : number of observations (n >= 2)

  The function Qn uses the procedures :
     whimed(a,pnW,n) : finds the weighted high median of an array
                      a of length n, using the array pnW (also of
                      length n) with positive integer weights.
     sort(x,n,y)    : sorts an array x of length n, and stores the
                      result in an array y (of size at least n)
     pull(a,n,k)    : finds the k-th order statistic of an array
                      a of length n

  When you want to use this function Qn in a program, you have to
  declare the next types in the mainprogram :
     type rvector = array [1..100] of real;
     type ivector = array [1..100] of integer;.
  You may change the number 100, which is the maximum size of n */

/*
void sort (double *pa, int n, double *pb)
{*/
/* Sorts an array a of length n and stores the result in the array b*/
/*
var jlv, jrv : ivector;
    i, jss, jndl, jr, jnc, j, jtwe : integer;
    xx, amm                        : real;
begin
  for i:= 1 to n do b[i] := a[i];
  jss := 1;
  jlv[1] := 1;
  jrv[1] := n;
  repeat
    jndl := jlv[jss];
    jr := jrv[jss];
    jss := jss - 1;
    repeat
      jnc := jndl;
      j := jr;
      jtwe := (jndl + jr) div 2;
      xx := b[jtwe];
      repeat
        while (b[jnc] < xx) do jnc := jnc + 1;
        while (xx < b[j]) do j := j -1;
        if (jnc <= j) then
          begin
            amm := b[jnc];
            b[jnc] := b[j];
            b[j] := amm;
            jnc := jnc + 1;
            j := j - 1
          end;
      until (jnc > j);
      if ((j - jndl) >= (jr - jnc)) then
        begin
          if (jndl < j) then
            begin
              jss := jss + 1;
              jlv[jss] := jndl;
              jrv[jss] := j
            end;
          jndl := jnc
        end
      else
        begin
          if (jnc < jr) then
            begin
              jss := jss + 1;
              jlv[jss] := jnc;
              jrv[jss] := jr
            end;
          jr := j
        end;
    until (jndl >= jr);
  until (jss = 0)
}*/


double pull(double *pdA, int n, int k)
{

// Finds the k-th order statistic of an array pdA of length n

	double *pdB = new double [n] ;				//	: rvector;
	double dAx, dBuffer ;		//: real;
	int l = 0, lr = n - 1, jnc, j ;

//	for i := 1 to n do
//		pdB[i] := pdA[i];
	memcpy (pdB, pdA, n * sizeof (double)) ;

	while (l < lr)
	{
		dAx = pdB[k];
		jnc = l;
		j = lr;
		while (jnc <= j)
		{
			while (pdB[jnc] < dAx)
				jnc++ ;
			while (pdB[j] > dAx)
				j-- ;

			if (jnc <= j)
			{
				dBuffer = pdB[jnc];
				pdB[jnc] = pdB[j];
				pdB[j] = dBuffer;
				jnc++ ;
				j-- ;
			}
		}
		if (j < k)
			l = jnc ;
		if (k < jnc)
			lr = j ;
	}

	dAx = pdB[k] ;		//	dAx wird ab hier nur mehr als temp  var verwendet..
	delete [] pdB ;

	return dAx ;

}


/* Algorithm to compute the weighted high median in O(n) time.

  The whimed is defined as the smallest a[j] such that the sum
  of the weights of all a[i] <= a[j] is strictly greater than
  half of the total weight.

  Parameters of this function :
     a : real array containing the observations
     n : number of observations
    pnW : array of integer weights of the observations.

  This function uses the function pull.

  The size of acand, iwcand must be at least n.*/

double whimed (double *pdA, int *pnW, int n)
{
	int i, nWtotal = 0, nWrest = 0, nWleft, nWmid, nWright, nKCand ;
    double dTrial ;
    double *pdAcand = new double [n] ;	//XX Length?!	//					: rvector;
	InitD (pdAcand, n, 0) ;
    int *pnWcand = new int [n] ;		//XX Length?!	//					: ivector;
	InitN (pnWcand, n, 0) ;
    bool bFound = false;				//					: boolean;

	double dRet = 0 ;

		//	for i = 1 to n do nWtotal = nWtotal + pnW[i];
	for (i = 0; i < n; i++)
		nWtotal += pnW[i];

	while (!bFound)
	{
			//	trial = pull(a,n,n div 2 + 1);
//		dTrial = pull(pdA,n,n / 2 + 1);
		dTrial = pull(pdA,n,n / 2);

		nWleft  = 0;
		nWmid = 0;
		nWright = 0;
			//	for i = 1 to n do
		for (i = 0; i < n; i++)
			if (pdA[i] < dTrial)
				nWleft += pnW[i] ;
			else if (pdA[i] > dTrial)
				nWright += pnW[i] ;
			else
				nWmid += pnW[i] ;
		
		if ((2 * nWrest + 2 * nWleft) > nWtotal)
		{
			nKCand = 0;
			for (i = 0; i < n; i++)
				if (pdA[i] < dTrial)
				{
					pdAcand[nKCand] = pdA[i] ;
					pnWcand[nKCand]= pnW[i] ;
					nKCand ++ ;
				}
			n = nKCand ;
		}
		else if ((2*nWrest+2*nWleft+2*nWmid) > nWtotal)
		{
			dRet = dTrial;
			bFound = true ;
		}
		else
		{
			nKCand = 0;
			for (i = 0; i < n; i++)
				if (pdA[i] > dTrial)
				{
					pdAcand[nKCand] = pdA[i];
					pnWcand[nKCand] = pnW[i] ;
					nKCand ++ ;
				}
			n = nKCand;
			nWrest += nWleft + nWmid ;
		}

		memcpy (pdA, pdAcand, n * sizeof (double)) ;
		memcpy (pnW, pnWcand, n * sizeof (int)) ;
/*		for (i = 0; i < n; i++)
		{
			pdA[i] = pdAcand[i] ;
			pnW[i] = pnWcand[i] ;
		}
*/
	}

    delete [] pdAcand ;
    delete [] pnWcand ;

	return dRet ;
}


/* Testsource

library (pcaPP)
r = matrix (rnorm (10000), 100, 100)
for (i in 1:1000)
  qna = qn(r)

*/

double myvec [] = {0.399, 0.994, 0.512, 0.844, 0.611, 0.857, 0.669, 0.872 } ;

	void Qn (double *pdX, int *pn, double *pdRet)
	{
		int n = *pn ;
		double *pdWork = new double [n]; //, *pdY = new double [n];							//XX Length?!	//	: rvector;

		InitD (pdWork, n, 0) ;
//		InitD (pdY, n, 0) ;
		
		int *pnLeft = new int [n],
			*pnRight = new int [n],
			*pnWeight = new int [n],
			*pnQ = new int [n],
			*pnP = new int [n];	//XX Length?!	//	: ivector;


		InitN (pnLeft, n, 0) ;
		InitN (pnRight, n, 0) ;
		InitN (pnWeight, n, 0) ;
		InitN (pnQ, n, 0) ;
		InitN (pnP, n, 0) ;

		double dn, dTrial, dQv = 0 ;										//	: real;
		int h, k, nKNew,nJHelp, nL, nR,nSumQ, nSumP, i, j, jj ;				//	: integer;
		bool bFound ;													//	: boolean;

			//h = n div 2 + 1;
		h = n / 2 + 1;
		k = h * (h - 1) / 2;

//		sort(pdX,n,pdY);
		double *pdY = new double [n] ;
		memcpy (pdY, pdX, n * sizeof (double) ) ;
		R_rsort (pdY, n) ;
//		doublesort (pdY, n) ;

			//	for (i = 0; i < n; i++)
		for (i = n-1; i >= 0; i--)
		{
			pnLeft[i] = n - i + 1;		//	 2 abziehen: 1 für den idx generell & 1 da ja mit geringeren indizes gerechnet wird..
			pnRight[i] = n ;
		}

		nJHelp = (n * (n + 1)) / 2;
		nKNew = k + nJHelp;
		nL = nJHelp;
		nR = n*n;
		bFound = false;
		while ((nR-nL > n) && (!bFound))
		{
			j = 0;
			for (i = 1; i < n; i++)
			{
				if (pnLeft[i] <= pnRight[i]) 
				{
					pnWeight[j] = pnRight[i] - pnLeft[i] + 1;
					nJHelp = pnLeft[i] + pnWeight[j] / 2;
					pdWork[j] = pdY[i] - pdY[n -  nJHelp];
					j++ ;
				}
			}
			dTrial = whimed(pdWork,pnWeight,j/* - 1*/);
			j = 0;
			for (i = n - 1; i >= 0; i--)
			{
				while ((j < n) && (pdY[i] - pdY[n - j - 1] < dTrial))
					j++ ;
				pnP[i] = j ;
			}
			j = n + 1;
			for (i = 0; i < n; i++)
			{
				while (pdY[i] - pdY[n - j + 1] > dTrial)
					j-- ;
				pnQ[i] = j ;
			}
			nSumP = 0;
			nSumQ = 0;

			//for (i = 0; i < n; i++)
			for (i = n-1; i >= 0; i--)
			{
				nSumP += pnP[i];
				nSumQ += pnQ[i] - 1 ;
			}
			if (nKNew <= nSumP) 
			{
				//	//	for i = 1 to n do
				//for (int i = n-1; i >= 0; i--)
				//	pnRight[i] = pnP[i] ;
				memcpy (pnRight, pnP, n* sizeof (int)) ;
				nR = nSumP ;
			}
			else if (nKNew > nSumQ) 
			{
				//	//for i = 1 to n
				//for (int i = n-1; i >= 0; i--)
				//	pnLeft[i] = pnQ[i];
				memcpy (pnLeft, pnQ, n* sizeof (int)) ;

				nL = nSumQ ;
			}
			else
			{
				dQv = dTrial;
				bFound = true ;
			}
		}
		if (!bFound) 
		{
			j = 0;
			for (i = 1; i < n; i++)
			{
				if (pnLeft[i] <= pnRight[i]) 
					for (jj = pnLeft[i]; jj <= pnRight[i]; jj++)
					{
						pdWork[j] = pdY[i] - pdY[n - jj];
						j++ ;
					}
			}
			dQv = pull(pdWork,j,nKNew-nL - 1) ;
		}

		if (n <=9)
			dn = myvec[n - 2] ;
		else if (n & 1)
			dn = n/(n+1.4) ;
		else
			dn = n/(n + 3.8) ;


/*		if (n <= 9) 
		{
			if (n=2)  dn = 0.399;
			if (n=3)  dn = 0.994;
			if (n=4)  dn = 0.512;
			if (n=5)  dn = 0.844;
			if (n=6)  dn = 0.611;
			if (n=7)  dn = 0.857;
			if (n=8)  dn = 0.669;
			if (n=9)  dn = 0.872
		}
		else 
		{
			if (n % 2 = 1)   dn = n/(n+1.4);
			if (n % 2 = 0)   dn = n/(n + 3.8)
		}*/

		delete [] pdY ;

		delete [] pdWork ;
		delete [] pnLeft ;
		delete [] pnRight ;
		delete [] pnWeight ;
		delete [] pnQ ;
		delete [] pnP ;

		*pdRet =  dn * 2.21914446598508
			* dQv ;
	}
