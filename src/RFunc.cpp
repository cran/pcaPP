/*
    IMat - "Intelligent" Matrix Classes
    Copyright (C) 2010 by Heinrich Fritz (heinrich_fritz@hotmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <math.h>
#include "RFunc.h"

//#include "IMat.h"


#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R.h>

	IMatD invert (IMatD A)
	{
		return invert_ref (A) ;
	}

	IMatD invert_ref (IMatD &A)
	{
		int p = A.ncol () ;
		IVec<int> v (p) ;
		IMatD B = IMatD::diag (p) ;

		int nInfo ;

		F77_CALL(dgesv) (&p, &p, A.FlatEdit (), &p, v.FlatEdit (), B.FlatEdit (), &p, &nInfo) ;

		return B ;
	}

	const IMatD &invert (const IMatD &A, const IMatD &B)
	{
		int p = A.ncol () ;
		IVec<int> v (p) ;
		B.setdiag () ;
		int nInfo ;
		F77_CALL(dgesv) (&p, &p, A.FlatEdit (), &p, v.FlatEdit (), B.FlatEdit (), &p, &nInfo) ;
		return B ;
	}

	BOOL eigen_sqr (const IMatD &A, const IVecD &vVal, const IMatD &mVec, BOOL bOrder)
	{
		int p = A.nrow () ;

		if ((DWORD) p != A.ncol ())
			return FALSE ;

		int nZ = 1, nw = p * 8, nInfo ; ;

		IVecD vvi (p) ;		//	imaginary part of eigenvalues
		IVecD vvw (nw) ;	//	work array

		IVecD cvVal (p) ;
		IMatD cmVec (p, p) ;

		IMatD temp (A, byVal);
//		temp <<= A ;

		F77_NAME(dgeev)("V", "N", &p, temp.FlatEdit (), &p, cvVal.FlatEdit(), vvi.FlatEdit (), cmVec.FlatEdit (), &p, NULL, &nZ, vvw.FlatEdit (), &nw, &nInfo) ;

		if (bOrder)
		{
			IVec <DWORD> ord (cvVal.order (TRUE)) ;
			vVal << cvVal (ISub(ord)) ;
			mVec << cmVec (ISub(), ISub(ord)) ;
		}
		else
		{
			vVal << cvVal ;
			mVec << cmVec ;
		}

		return TRUE ;
	}


	BOOL eigen_val (const IMatD &A, IVecD &vVal, IMatD &mTemp, IVecD &vTemp)
	{
		int p = A.ncol (), n = A.nrow () ;

		mTemp.Reshape (n, p) << A ;
		vVal.Reshape (min (p, n)) ;

		int nInfo ;
		int nWork = -1 ;
		int nZ = 1 ;
		double dWork ;

		F77_NAME(dgesvd)("N", "N", &n, &p, NULL, &n, NULL, &dWork, &nZ, &dWork, &nZ, &dWork, &nWork, &nInfo) ;

		vTemp.Reshape ((DWORD) (nWork = int (ceil (dWork)))) ;

		F77_NAME(dgesvd)("N", "N", &n, &p, mTemp.FlatEdit (), &n, vVal.FlatEdit (), &dWork, &nZ, &dWork, &nZ, vTemp.FlatEdit (), &nWork, &nInfo) ;

		return TRUE ;
	}



	BOOL eigen_val_sqr (const IMatD &A, const IVecD &vVal, BOOL bOrder)
	{
		int p = A.nrow () ;

		if ((DWORD) p != A.ncol ())
			return FALSE ;

		int nZ = 1, nw = p * 8, nInfo ; ;

		IVecD vvi (p) ;		//	imaginary part of eigenvalues
		IVecD vvw (nw) ;	//	work array

		IVecD cvVal (p) ;

		IMatD temp (A, byVal);

		F77_NAME(dgeev)("N", "N", &p, temp.FlatEdit (), &p, cvVal.FlatEdit (), vvi.FlatEdit (), NULL, &nZ, NULL, &nZ, vvw.FlatEdit (), &nw, &nInfo) ;

		IVec <DWORD> ord (cvVal.order (TRUE)) ;

		if (bOrder)
			vVal << cvVal (ISub(ord)) ;
		else
			vVal << cvVal ;

		return TRUE ;
	}

	void prod (const IVecD &vec, double &val)
	{
		DWORD v ;
		val = 1 ;

		for (v = vec.size () - 1; v != (DWORD) -1; v--)
			val *= vec(v) ;
	}

	void sum (const IVecD &vec, double &val)
	{
		DWORD v ;
		val = 0 ;

		for (v = vec.size () - 1; v != (DWORD) -1; v--)
			val += vec(v) ;
	}

	void sumsq (const IVecD &vec, double &val)
	{
		DWORD v ;
		val = 0 ;

		
		for (v = vec.size () - 1; v != (DWORD) -1; v--)
		{
			double &dCur = vec(v) ;
			val += dCur * dCur ;
		}
	}

	void sqrtsumsq (const IVecD &vec, double &val)
	{
		DWORD v ;
		val = 0 ;

		
		for (v = vec.size () - 1; v != (DWORD) -1; v--)
		{
			double &dCur = vec(v) ;
			val += dCur * dCur ;
		}
		val = sqrt (val) ;
	}

	void sumabs (const IVecD &vec, double &val)
	{
		DWORD v ;
		val = 0 ;
		for (v = vec.size () - 1; v != (DWORD) -1; v--)
			val += fabs (vec(v)) ;
	}

//	void min (const IVecD &vec, double &val)
	void min (const IVecConstD &vec, double &val)
	{
		val = vec(0) ;
		DWORD v = vec.size () - 1 ;
		for (; v; v--)
		{
			const double d = vec(v) ;
			if (val > d)
				val = d ;
		}
	}

	void max (const IVecConstD &vec, double &val)
	{
		val = vec(0) ;
		DWORD v = vec.size () - 1 ;
		for (; v; v--)
		{
			const double d = vec(v) ;
			if (val < d)
				val = d ;
		}
	}

	double &maxDiag (const IMatD &mat, double &dMaxDiag)
	{
		if (!mat.ncol () || !mat.nrow ())	//	matrix has no extend
			return dMaxDiag = R_NaN ;		

		dMaxDiag = mat(0, 0) ;

		DWORD i ;
		for (i = min (mat.ncol (), mat.nrow ()) - 1; i != 0; i--)
			setmax (mat(i, i), dMaxDiag) ;
		return dMaxDiag ;
	}

	void Pmax (const IVecConstD &vec, double &val, DWORD dwEnd)
	{
		val = vec(0) ;
		DWORD v = dwEnd - 1 ;
		for (; v; v--)
			setmax (vec(v), val) ;
	}

	void iter_row (const IMatD &mat, const IVecD &vec, void (*func) (const IVecD &, double &))
	{
		ASSERT (vec.size () >= mat.nrow ()) ;
		DWORD i ;

		for (i = mat.nrow () - 1; i != (DWORD) -1; i--)
			func (mat.GetRow (i), vec(i)) ;
	}

	void iter_col (const IMatD &mat, const IVecD &vec, void (*func) (const IVecD &, double &))
	{
		ASSERT (vec.size () >= mat.ncol ()) ;
		DWORD i ;
		for (i = mat.ncol () - 1; i != (DWORD) -1; i--)
		{
			IVecD col = mat.GetCol (i) ;
			func (col, vec(i)) ;
		}
	}

	void cov (const IMatD &mat, const IMatD &ret)
	{
		IMatD cmat (mat.byrow () -  colMeans (mat)) ;
		matmultmat (t(cmat), cmat, ret) ;
		ret *= 1 / ((double) mat.nrow () - 1) ;
//		ret = (t(cmat) ->* cmat) / (mat.nrow () - 1) ;
	}

/* Equal probability sampling; without-replacement case */

	DWORD *SampleNoReplace(int k, int n, DWORD *y, DWORD *x)
	{
		int i, j;
		for (i = n - 1; i != -1; i--)
			x[i] = i ;

		for (i = 0; i < k; i++)
//		for (i = k - 1; i != -1; i--)
		{
			j = int (n * unif_rand()) ;
			y[i] = x[j] ;// + 1;
			x[j] = x[--n];
		}
		return x ;
	}

	void SampleNoReplace(int k, int n, const IVec<DWORD> &vec)
	{
		delete [] SampleNoReplace (k, n, vec.FlatEdit (), new DWORD [n]) ;
	}

	IVec<DWORD> SampleNoReplace (int k, int n)
	{
		IVec<DWORD> v (k) ;
		SampleNoReplace (k, n, v) ;
		return v ;
	}

	void rank (IVecD &v, IVecDW &r)
	{
		r.Reshape(v.size ()) ;
		rsort_with_index (v.FlatEdit (), (int *) (DWORD *) r.FlatEdit (), v.size ()) ;
	}

	const IVecD &sort (const IVecD &v, BOOL bDesc)
	{
		R_rsort (v.FlatEdit (), v.size ()) ;
		return v ;

//		if (bDesc)
//			sorted.InvertOrder () ;
	}

	void Rprin1212tVec (const IVecD &vec)
	{
		DWORD i ;
		for (i = 0; i < vec.size (); i++)
			Rprintf ("%.5f\r\n", vec (i)) ;
	}


	void RtprintVec (const IVecD &vec, LPCTSTR szText, BOOL bAppendCRLF)
	{
		if (szText)
			Rprintf (szText) ;

		DWORD i ;
		for (i = 0; i < vec.size (); i++)
			Rprintf ("%.10f\t", (double) vec (i)) ;
		if (bAppendCRLF)
			Rprintf ("\r\n") ;
	}

	void RtprintVec (const IVecDW &vec, LPCTSTR szText, BOOL bAppendCRLF)
	{
		if (szText)
			Rprintf (szText) ;
		DWORD i ;
		for (i = 0; i < vec.size (); i++)
			Rprintf ("%d ", vec (i)) ;
		if (bAppendCRLF)
			Rprintf ("\r\n") ;
	}

	void RprintMat (const IMatD &mat, LPCTSTR szText)
	{
		if (szText)
			Rprintf (szText) ;
		DWORD i, j ;
		for (i = 0; i < mat.nrow (); i++)
		{
			for (j = 0; j < mat.ncol (); j++)
				Rprintf ("%.10e \t", mat (i, j)) ;
			Rprintf ("\r\n") ;
		}
	}

//static void rPsort2(double *x, int lo, int hi, int k) ;

	double median (const IVecD &v)
	{
		IVecD temp (v.size ()) ;
		temp << v ;
		return median_raw (temp) ;
	}

	double median (const IVecD &v, IVecD &temp)
	{
		temp << v ;
		return median_raw (temp) ;
	}

	double median_raw (IVecD &v)
	{
		DWORD dwHalf = v.size () >> 1 ;
		rPsort (v.FlatEdit (), v.size(), dwHalf) ;
		if (v.size () % 2)	//	odd -> just get 1 value
			return v(dwHalf) ;
		return (Pmax (v, dwHalf) + v(dwHalf)) / 2 ;
	}

	void ColSumWeighted (const IMatD &m, const IVecD &v, IVecD &s)
	{
		ASSERT (m.nrow () == v.size ()) ;

		s.Reshape (m.ncol ()) ;
		s.Reset (0) ;

		DWORD r, c, rm1 = m.nrow () - 1 ;

		for (c = m.ncol () - 1 ;c != (DWORD) -1; c--)
			for (r = rm1; r != (DWORD) -1; r--)
				s (c) += v(r) * m(r, c) ;
	}

	DWORD which_max_abs (const IVecD &v)
	{
		DWORD dwIdx = 0 ;
		double dMax = fabs (v(0)) ;
		DWORD i ;
		
		for (i = v.size () - 1; i; i--)
			if (setmax (fabs (v(i)), dMax))
				dwIdx = i  ;

		return i ;
	}

	DWORD first_idx_NZ (const IVecD &v, double dZeroTol)
	{
		DWORD i ;
		for (i = v.size () - 1; i != (DWORD) -1; i--)
			if (fabs (v(i)) > dZeroTol)
				return i ;
		return -1 ;
	}
