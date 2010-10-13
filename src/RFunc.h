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


#ifdef _MSC_VER
	#define EXPORT extern "C" __declspec(dllexport)
#else
	#define EXPORT extern "C"
#endif //#ifdef WIN32


#include "IMat.h"



	IMatD invert_ref (IMatD &A) ;
	BOOL eigen_sqr (const IMatD &A, const IVecD &vVal, const IMatD &mVec, BOOL bOrder = TRUE) ;
	BOOL eigen_val_sqr (const IMatD &A, const IVecD &vVal, BOOL bOrder = TRUE) ;
	BOOL eigen_val (const IMatD &A, IVecD &vVal, IMatD &mTemp, IVecD &vTemp) ;

	IMatD invert (IMatD A) ;
	const IMatD &invert (const IMatD &A, const IMatD &B) ;

	void sumsq (const IVecD &vec, double &val) ;
	inline double sumsq (const IVecD &vec) { double d ; sumsq (vec, d) ; return d ; }

	void sumabs (const IVecD &vec, double &val) ;
	inline double sumabs (const IVecD &vec) { double d ; sumabs (vec, d) ; return d ; }

	void sqrtsumsq (const IVecD &vec, double &val) ;
	inline double sqrtsumsq (const IVecD &vec) { double d ; sqrtsumsq (vec, d) ; return d ; }

	void sum (const IVecD &vec, double &val) ;
	inline double sum (const IVecD &vec) { double d ; sum (vec, d) ; return d ; }

	void prod (const IVecD &vec, double &val) ;
	inline double prod (const IVecD &vec) { double d ; prod (vec, d) ; return d ; }

	inline void mean (const IVecD &vec, double &val)	{ val = sum (vec) / vec.size () ; }
	inline double mean (const IVecD &vec) { double d ; mean (vec, d) ; return d ; }

	void min (const IVecConstD &vec, double &val) ;
	inline double min (const IVecConstD &vec) { double d ; min (vec, d) ; return d ; }

	void max (const IVecConstD &vec, double &val) ;
	inline double max (const IVecConstD &vec) { double d ; max (vec, d) ; return d ; }

	void inline min_ (const IVecD &vec, double &val) { min (IVecConstD (vec), val) ; }
	void inline max_ (const IVecD &vec, double &val) { max (IVecConstD (vec), val) ; }

	void Pmax (const IVecConstD &vec, double &val, DWORD dwEnd) ;
	inline double Pmax (const IVecConstD &vec, DWORD dwEnd) { double d ; Pmax (vec, d, dwEnd) ; return d ; }


	template <class T> inline bool setmax (const T &val, T &max)
	{
		if (val < max)
			return false ;
		max = val ;
		return true ;
	}

	template <class T> inline bool setmin (const T &val, T &min)
	{
		if (val > min)
			return false ;
		min = val ;
		return true ;
	}

	template <class T> inline void setminmax (const T &val, T &min, T&max) 
	{
		if (val < min)	min = val ;
		else if (val > max)	max = val ;
	}

/*	template <class T> DWORD count (const ITens_ <T> &t, DWORD &dwCount) ;
	{
		dwCount = 0 ;
		ITens_<T>::t_Iter it = t.dim () - 1 ;

		do
		{
			if (t(it) != 0)
				dwCount++ ;
		} while (t.It_mm (it)) ;

		return dwCount ;
	}
*/

	void iter_row (const IMatD &mat, const IVecD &vec, void (*func) (const IVecD &, double &)) ;
	void iter_col (const IMatD &mat, const IVecD &vec, void (*func) (const IVecD &, double &)) ;

	
	double &maxDiag (const IMatD &mat, double &dMaxDiag) ;
	inline double maxDiag (const IMatD &mat) { double dMaxDiag ; return maxDiag (mat, dMaxDiag) ;}
	

	void cov (const IMatD &mat, const IMat<double> &ret) ;

	inline IMatD cov (const IMatD &mat)
	{
		IMat<double> ret (mat.ncol (), mat.ncol ()) ;
		cov (mat, ret) ;
		return ret ;
	}

	inline void colSums (const IMatD &A, const IVecD &vec)	{ iter_col (A, vec, sum) ; }
	inline void rowSums (const IMatD &A, const IVecD &vec)	{ iter_row (A, vec, sum) ; }
	inline  IVecD colSums (const IMatD &A) { IVecD v (A.ncol ()) ; colSums (A, v) ; return v ; }
	inline IVecD rowSums (const IMatD &A) { IVecD v (A.nrow ()) ; rowSums (A, v) ; return v ; }

	inline void colProds (const IMatD &A, const IVecD &vec)	{ iter_col (A, vec, prod) ; }
	inline void rowProds (const IMatD &A, const IVecD &vec)	{ iter_row (A, vec, prod) ; }
	inline  IVecD colProds (const IMatD &A) { IVecD v (A.ncol ()) ; colProds (A, v) ; return v ; }
	inline IVecD rowProds (const IMatD &A) { IVecD v (A.nrow ()) ; rowProds (A, v) ; return v ; }

	inline void colMins (const IMatD &A, const IVecD &vec)	{ iter_col (A, vec, min_) ; }
	inline void rowMins (const IMatD &A, const IVecD &vec)	{ iter_row (A, vec, min_) ; }
	inline  IVecD colMins (const IMatD &A) { IVecD v (A.ncol ()) ; colMins (A, v) ; return v ; }
	inline IVecD rowMins (const IMatD &A) { IVecD v (A.nrow ()) ; rowMins (A, v) ; return v ; }

	inline void SqrtrowSumSqs (const IMatD &A, const IVecD &vec)	{ iter_row (A, vec, sqrtsumsq) ; }

	void ColSumWeighted (const IMatD &m, const IVecD &v, IVecD &s) ;

	DWORD which_max_abs (const IVecD &v) ;
	DWORD first_idx_NZ (const IVecD &v, double dZeroTol = 1e-16) ;

	inline void colMeans (const IMatD &A, const IVecD &vec)	{ iter_col (A, vec, mean) ; }
	inline void rowMeans (const IMatD &A, const IVecD &vec)	{ iter_row (A, vec, mean) ; }
	inline  IVecD colMeans (const IMatD &A) { IVecD v (A.ncol ()) ; colMeans (A, v) ; return v ; }
	inline IVecD rowMeans (const IMatD &A) { IVecD v (A.nrow ()) ; rowMeans (A, v) ; return v ; }


	DWORD *SampleNoReplace(int k, int n, DWORD *y, DWORD *x) ;
	void SampleNoReplace(int k, int n, const IVec<DWORD> &vec) ;
	IVec<DWORD> SampleNoReplace (int k, int n) ;
	const IVecD &sort (const IVecD &v, BOOL bDesc = FALSE) ;
	void rank (IVecD &v, IVecDW &r) ;
/*
	template <class T, class D, class U>
	BOOL Equal (const ITens <T, D> t1, const ITens <U, D> t2, BOOL bEqualSize = FALSE)
	{
		if (!bEqualSize)
		{
			DWORD d ;
			for (d = t1.ndim () - 1; d != (DWORD) -1; d--)
				if (t1.size (d) != t2.size (d))
					return FALSE ;
		}
	}
*/

		void RprintVec (const IVecD &vec) ;
		void RtprintVec (const IVecD &vec, LPCTSTR szText = NULL, BOOL bAppendCRLF = TRUE) ;
		void RtprintVec (const IVecDW &vec, LPCTSTR szText = NULL, BOOL bAppendCRLF = TRUE) ;
		void RprintMat (const IMatD &mat, LPCTSTR szText = NULL) ;

		double median (const IVecD &v) ;
		double median (const IVecD &v, IVecD &temp) ;
		double median_raw (IVecD &v) ;


