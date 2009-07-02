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


#pragma once

#include "ITens.h"
#include <stdio.h>

template <class T> class IMatByRow ;


	template <class TA, class TB, class T>
	BOOL matmultmat (const TA &a, const TB &b, const IMat<T> &res) ;

	template <class TA, class TB, class TTC>
	BOOL vecmultmat (TA &a, TB &b, IMat<TTC> &res) ;

	template <class TA, class TB, class TTC>
	BOOL matmultvec (TA &a, TB &b, IMat<TTC> &res) ;

	template <class T>
	class OrderContainer
	{
		public:

		static inline int compareDesc (OrderContainer<T> *p1, OrderContainer<T> *p2) { return -compare (p1, p2) ;}

		static int compare (OrderContainer<T> *p1, OrderContainer<T> *p2)
		{
			if (p1->m_val < p2->m_val)
				return -1 ;
			if (p1->m_val > p2->m_val)
				return 1 ;
			return 0 ;
		}

		T m_val ;
		DWORD m_dwIdx ;
	} ;

	template <class T>
	class SortFunc
	{
		public:

		static inline int compareDesc (T *p1, T *p2) { return -compare (p1, p2) ;}

		static int compare (T *p1, T *p2)
		{
			if (*p1 < *p2)
				return -1 ;
			if (*p1 > *p2)
				return 1 ;
			return 0 ;
		}

	} ;

///////////////////////////////////////////////////////////////////////////////
//	IVecConst
///////////////////////////////////////////////////////////////////////////////

	template <class T>
	class IVecConst : public ITensConst <T, IDimsVec <T> >
	{
		protected:
			typedef ITensConst <T, IDimsVec <T> > t_base ;
		public:

		inline T & operator ()	(const DWORD &c)	const { return t_base::m_pDataRef->GetData () [t_base::GetIdx (c)] ; }
	} ;

///////////////////////////////////////////////////////////////////////////////
//	IVecEdit
///////////////////////////////////////////////////////////////////////////////

/*	template <class T>
	class IVecEdit : public ITensEdit_ <T, IDimsVec <T> >
	{
	protected:
		typedef CDataRef <T> m_reftype ;

		typedef IDimsVec <T> t_dims ;
		typedef ITensEdit_ <T, IDimsVec <T> > t_tensedit, t_base ;
	public:

		
//		inline IVecEdit <T> &operator = (const ITensEdit_<T, IDimsVec <T> > &tensedit) const	{ return FC_As::OpTT (*this, tensedit) ; }

		IVecEdit () {} ;
		IVecEdit (const IVecEdit <T> &vecedit) : t_tensedit (vecedit ) {}
		IVecEdit (m_reftype *&pDataRef, DWORD dwOffset, IDim *pDim)	: t_tensedit (pDataRef, dwOffset, pDim) {}
		IVecEdit (m_reftype *&pDataRef, const t_dims &ida)	: t_tensedit (pDataRef, ida) {}

		inline T & operator ()	(DWORD c) const { return Element (c) ; }
		inline T & Element	(DWORD c) const { return t_base::m_pDataRef->GetData () [t_dims::GetIdx (c)] ; }

		inline IVec<T> operator ()	(const ISub_ &v) const { return IVecEdit<T> (t_base::m_pDataRef, t_base::m_dwOffset, IDim::Create (v, t_base::m_apDim[0])) ; }

//		inline const IVecEdit <T> &operator = (const IVecEdit <T> &vec)	const { return FC_As::OpVV (*this, vec) ; }
//		inline const IVecEdit <T> &operator = (const IVec <T> &vec)		const { FC_As::OpVV (*this, vec) ; return *this ;}
//		inline const IVecEdit <T> &operator = (const ITens (T, 1) &tens)	const { FC_As::OpTT (*this, tens) ; return *this ;}

		operator const IVecConst <T> & () const { return * (IVecConst <T> *) this ; }
	} ;
*/
///////////////////////////////////////////////////////////////////////////////
//	IVec
///////////////////////////////////////////////////////////////////////////////

	template <class T>
	class IVec : public ITens_ <T, IDimsVec <T> >
	{
//		friend class IVecEdit <T> ;

		typedef ITens_ <T, IDimsVec <T> > t_tens, t_base ;
		typedef IDimsVec <T> t_dims ;

	public:
//		typedef IVecEdit <T> m_edittype ;
		typedef CDataRef <T> m_reftype ;
		typedef T m_datatype ;

		IVec () {}
		IVec (DWORD dwOffset, IDim *pDim, CDataRef <T> *pDataRef) : ITens_ <T, t_dims> (dwOffset, pDataRef, pDim/*, NULL*/)	{ }

		IVec (DWORD dwSize, T *pData = NULL, SetDataFlag sdf = sdf_modeReference)  : ITens_ <T, t_dims> (pData, sdf, dwSize/*, -1*/) {}
		IVec (const IVec<T> &vec, const createMode cm = byRef) : ITens_<T, t_dims> (vec, cm) {}

		IVec (m_reftype *pref, DWORD dwOffset, IDim *pdim) : ITens_ <T, t_dims> (dwOffset, pref, pdim/*, NULL*/) {}
		IVec (CDataRef <T> *pDataRef, const IDimsVec<T> &ida) : t_base (pDataRef, ida) {}

		inline DWORD size () const { return t_base::m_apDim[0]->GetSize ();}

		inline T &operator () (DWORD c) const { return t_base::m_pDataRef->GetData () [t_base::GetIdx (c)] ; }
		inline IVec<T> operator () (const ISub_ &v)	const { return IVec<T> (t_base::m_pDataRef, t_base::m_dwOffset, IDim::Create (v, t_base::m_apDim[0])) ; }

		void Create (const DWORD &c, T *pData = NULL, SetDataFlag sdf = sdf_modeReference) 
		{ 
			t_base::m_pDataRef->ReCreate (t_dims::ReCreate (c), t_base::m_pDataRef, pData, sdf, this) ;
		}

//		inline operator IVecEdit<T> () { return t_tens::Edit () ; }
//		inline IVecEdit <T> Edit () { return t_base::Edit () ; }

		void print ()
		{
			DWORD i ;
			for (i = 0; i < size (); i++)
				::printf ("%.0f\r\n", (double) t_tens::Element (i)) ;
		}

		void tprint ()
		{
			DWORD i ;
			for (i = 0; i < size (); i++)
				::printf ("%.0f\t", (double) t_tens::Element (i)) ;
			::printf ("\r\n") ;
		}

		IVec <DWORD> order (BOOL bDecreasing = FALSE)
		{
			IVec <DWORD> ret (size ()) ;
			order (ret, bDecreasing) ;
			return ret ;
		}

		void order (const IVec<DWORD> &ord, BOOL bDecreasing = FALSE)
		{
			OrderContainer<T> *pOrd = new OrderContainer<T> [size ()] ;

			DWORD v ;
			for (v = size () - 1; v != (DWORD) -1; v--)
			{
				pOrd[v].m_dwIdx = v ;
				pOrd[v].m_val = Element (v) ;
			}

			if (bDecreasing)
				qsort (pOrd, size (), sizeof (OrderContainer<T>), (int (*)(const void*, const void*)) OrderContainer<T>::compareDesc) ;
			else
				qsort (pOrd, size (), sizeof (OrderContainer<T>), (int (*)(const void*, const void*)) OrderContainer<T>::compare) ;


			for (v = size () - 1; v != (DWORD) -1; v--)
				ord(v) = pOrd[v].m_dwIdx ;

			delete [] pOrd ;
		}

		T *Detach (DWORD &dwSize = *(DWORD *) NULL)
		{
			if (!t_base::IsFlat ())
				return NULL ;

			if (&dwSize)
				dwSize = size () ;
			IDim::Empty (t_base::m_apDim[0]) ;

			return ::Detach (t_base::m_pDataRef) ;
		}

		IMPL_OPERATOR_VEC(+, += , FC::FC_plus, IVec, T)
		IMPL_OPERATOR_VEC(-, -= , FC::FC_minus, IVec, T)
		IMPL_OPERATOR_VEC(/, /= , FC::FC_divide, IVec, T)
		IMPL_OPERATOR_VEC(*, *= , FC::FC_multiply, IVec, T)
		IMPL_OPERATOR_VEC(^, ^= , FC::FC_pow, IVec, T)

		IMPL_OPERATOR_VEC_NOAS (>, FC::FC_greater, CVecL, bool)
		IMPL_OPERATOR_VEC_NOAS (>=, FC::FC_greater_equal, CVecL, bool)
		IMPL_OPERATOR_VEC_NOAS (==, FC::FC_equal, CVecL, bool)
		IMPL_OPERATOR_VEC_NOAS (<, FC::FC_less, CVecL, bool)
		IMPL_OPERATOR_VEC_NOAS (<=, FC::FC_less_equal, CVecL, bool)

		inline IVec<bool> operator - () const { return FC_ElOp<FC::FC_minus_u, bool>::OpV (*this); }

		operator const IVecConst <T> & () const { return * (IVecConst <T> *) this ; }

		IVec <T> Reshape (const DWORD r)
		{
			if (r > t_base::m_pDataRef->GetSize ())
				return Recreate (r)  ;

			IDim::ReCreate (r, 1, t_dims::m_apDim[0]) ;

			return *this ;
		}

		IVec <T> &Recreate (DWORD c)
		{
//			IDim::ReCreate (c, r, t_dims::m_apDim[0]) ;

			t_base::Recreate (IIter <1> (c)) ;
			return *this ;
		}

	protected:
		inline T & Element (DWORD c) const { return t_base::m_pDataRef->GetItem ( t_base::GetIdx (c)) ; }	//	this must only be called, when in edit (or single owner) mode		
	} ;

	typedef IVec<double> IVecD ;
//	typedef IVecEdit<double> IVecEditD ;
	typedef IVecConst<double> IVecConstD ;
	typedef CSimpleArray <IVec <double> > IVecArrayD  ;
//	typedef CSimpleArray <IVecEdit <double> > IVecEditArrayD  ;

	typedef IVec<DWORD> IVecDW ;
	typedef IVecConst<DWORD> IVecConstDW ;
	typedef CSimpleArray <IVec <DWORD> > IVecArrayDW  ;
	typedef IVec<int> IVecN ;

	typedef IVec<bool> IVecB ;

///////////////////////////////////////////////////////////////////////////////
//	IVecL
///////////////////////////////////////////////////////////////////////////////

	template <class T>
	class IVecL : public IVec <T>
	{
		public:
		IMPL_OPERATOR_VEC(|, |= , FC::FC_or,	IVecL, T)
		IMPL_OPERATOR_VEC(&, &= , FC::FC_and,	IVecL, T)
		IMPL_OPERATOR_VEC_NOAS(||, FC::FC_OR,	IVecL, bool)
		IMPL_OPERATOR_VEC_NOAS(&&, FC::FC_AND,	IVecL, bool)

		inline IVecL<bool> operator ! () const { return FC_ElOp<FC::FC_NOT, bool>::OpV (*this); }
		inline IVecL<T> operator ~ () const { return FC_ElOp<FC::FC_not, bool>::OpV (*this); }
	} ;

///////////////////////////////////////////////////////////////////////////////
//	IMatConst
///////////////////////////////////////////////////////////////////////////////

	template <class T>
	class IMatConst : public ITensConst <T, IDimsMat <T> >
	{
		protected:
			typedef ITensConst <T, IDimsMat <T> > t_base ;
		public:

		inline T & operator ()	(const DWORD &r, const DWORD &c)	const { return t_base::m_pDataRef->GetData () [t_base::GetIdx (r, c)] ; }
	} ;

///////////////////////////////////////////////////////////////////////////////
//	IMatEdit
///////////////////////////////////////////////////////////////////////////////

/*	template <class T>
	class IMatEdit : public ITensEdit_ <T, IDimsMat <T> > //CEditBase <CDataRef <T> >
	{
		typedef IDimsMat <T> t_dims ;
		typedef ITensEdit_ <T, t_dims > t_tensedit, t_base ;
		
		
	public:
		typedef CDataRef <T> m_reftype ;
		typedef T m_datatype ;
		typedef IMatEdit <T> m_edittype, t_editbase ;

//		IMatEdit (const IMatEdit <T> &imatedit) : CEditBase <m_reftype> (imatedit), m_pDimR (NULL), m_pDimC (NULL), m_pdData (imatedit.m_pdData)
//		{
//			imatedit.m_pDimR->Reference (m_pDimR) ;
//			imatedit.m_pDimC->Reference (m_pDimC) ;
//		}
//
//		IMatEdit (m_reftype *&pDataRef, IDim *pDimR, IDim *pDimC) : CEditBase <m_reftype> (pDataRef), m_pDimR (NULL), m_pDimC (NULL), m_pdData (pDataRef->GetData ())
//		{
//			pDimR->Reference (m_pDimR) ;
//			pDimC->Reference (m_pDimC) ;
//		}


		IMatEdit (const IMatEdit <T> &vecedit) : t_tensedit (vecedit ) {}
		IMatEdit (m_reftype *&pDataRef, DWORD dwOffset, IDim *pDim)	: t_tensedit (pDataRef, dwOffset, pDim) {}
		IMatEdit (m_reftype *&pDataRef, const t_dims &ida)	: t_tensedit (pDataRef, ida) {}


		inline T & operator ()	(DWORD r, DWORD c) const { return Element (r, c) ; }
		inline T & Element		(DWORD r, DWORD c) const { return t_editbase::m_pDataRef->GetData () [t_base::GetIdx (r, c)] ; }

		

//		inline DWORD nrow () { return t_base::m_apDim[0]->GetSize () ; }
//		inline DWORD ncol () { return t_base::m_apDim[1]->GetSize () ; }

		inline operator T * () { return t_tensedit::Data () ; }

		inline IMatEdit<T> &Edit () { return *this ; }

		void Reset (T val = 0)
		{
			DWORD i, j ;
			for (i = t_base::nrow () - 1; i != (DWORD) -1; i--)
				for (j = t_base::ncol () - 1; j != (DWORD) -1; j--)
					Element (i, j) = val ;
		}

//		inline const IMatEdit <T> &operator = (const IMatEdit <T> &mat) const	{ return FC_As::OpMM (*this, mat) ; }
//		inline const IMatEdit <T> &operator = (const IMat <T> &mat) const { return FC_As::OpMM (*this, mat) ; }
//		template <class U> inline IMatEdit <T> & operator = (const IMat<U> mat) { return FC_As::OpMM (*this, mat) ; }

		operator const IMatConst <T> & () const { return * (IMatConst <T> *) this ; }

	protected:

	} ;
*/
///////////////////////////////////////////////////////////////////////////////
//	IMat
///////////////////////////////////////////////////////////////////////////////
/*inline IDim *PrintDim (IDim *pDim, const char *sz)
{
		Rprintf ("handling dim %08p of size %d", pDim, pDim->GetSize ()) ;

		if (sz) Rprintf (sz) ;
		Rprintf ("\r\n") ;
		return pDim ;

}*/
	class TClust ;
	template <class T>
	class IMat : public ITens_ <T, IDimsMat <T > >
	{
		friend class TClust ;
		friend class TClust1 ;
		typedef IDimsMat <T> t_dims ;
		typedef ITens_ <T, t_dims > t_tens, t_base ;
		typedef typename t_base::t_Iter t_Iter ;

	public:
		typedef CDataRef <T> m_reftype ;
		typedef T m_datatype ;
		typedef ITensFlatEdit < T, IDimsMat <T > > t_flatedit	;

		IMat () {} 

		IMat (m_reftype *pref, DWORD dwOffset, IDim *pdimr, IDim *pdimc) : t_tens  (dwOffset, pref, pdimr, pdimc) { }

		IMat (DWORD dwNR, DWORD dwNC, T *pData = NULL, BOOL bCM = TRUE, SetDataFlag sdf = sdf_modeReference) : t_tens (pData, sdf, bCM ? dwNR : dwNC, bCM ? dwNC : dwNR)
		{
			if (!bCM)
				t () ;
		}

		void Create (const DWORD &r, const DWORD &c, T *pData = NULL, SetDataFlag sdf = sdf_modeReference) 
		{ 
			t_base::m_pDataRef->ReCreate (t_dims::ReCreate (t_Iter (r, c)), t_base::m_pDataRef, pData, sdf, this) ;
		}


		IMat (const IMat<T> &mat, const createMode cm = byRef) : t_tens (mat, cm) {}

		IMat (CDataRef <T> *pDataRef, const IDimsMat<T> &ida) : t_base (pDataRef, ida) {}

		inline IVec<T> GetCol (DWORD dwIdx) const { return IVec<T>(t_tens::m_dwOffset + t_tens::m_apDim[1]->Get_NC (dwIdx), t_tens::m_apDim[0], t_tens::m_pDataRef) ; }
		inline IVec<T> GetRow (DWORD dwIdx) const { return IVec<T> (t_tens::m_dwOffset + t_tens::m_apDim[0]->Get_NC (dwIdx), t_tens::m_apDim[1], t_tens::m_pDataRef) ; }

		m_reftype *GetReference (m_reftype *&pref) const { return t_base::m_pDataRef->Reference (pref) ; }

	//	mathematical operations
		IMat<T> operator () (const ISub_ &r,  const ISub_ &c) const { return IMat<T> (t_tens::m_pDataRef, t_base::m_dwOffset, IDim::Create (r, t_base::m_apDim[0]), IDim::Create (c, t_base::m_apDim[1])) ; }

		inline T &operator () (DWORD r, DWORD c) const { return Element (r, c) ; }

		IMat &t ()
		{
//			ASSERT (!IsLocked ()) ;
			swap (&t_tens::Dim (0), &t_tens::Dim (1)) ;
			return *this ;
		}

		template <class U>
		inline IMat <T> operator ->* (IVec<U> &vec) const 
		{
			IMat <T> res ;
			matmultvec (*this, vec, res) ;
			return res ;
		}

		template <class U>
		inline IMat <T> operator ->* (const IMat<U> &mat) const
		{
			IMat <T> res (t_base::nrow (), mat.ncol ()) ;

			matmultmat (*this, mat, res) ;
			return res ;
		}

		

		IMPL_OPERATOR_MAT (+, +=, FC::FC_plus,		IMat, T) ;
		IMPL_OPERATOR_MAT (-, -=, FC::FC_minus,		IMat, T) ;
		IMPL_OPERATOR_MAT (/, /=, FC::FC_divide,	IMat, T) ;
		IMPL_OPERATOR_MAT (*, *=, FC::FC_multiply,	IMat, T) ;
		IMPL_OPERATOR_MAT (^, ^=, FC::FC_pow,		IMat, T) ;

		IMPL_OPERATOR_MAT_NOAS (>, FC::FC_greater, CMatL, bool)
		IMPL_OPERATOR_MAT_NOAS (>=, FC::FC_greater_equal, CMatL, bool)
		IMPL_OPERATOR_MAT_NOAS (==, FC::FC_equal, CMatL, bool)
		IMPL_OPERATOR_MAT_NOAS (<, FC::FC_less, CMatL, bool)
		IMPL_OPERATOR_MAT_NOAS (<=, FC::FC_less_equal, CMatL, bool)

		inline IMat<bool> operator - () const { return FC_ElOp<FC::FC_minus_u, bool>::OpV (*this); }

	//	Operational Items

		void print () 
		{
			DWORD i, j ;
			for (i = 0; i < t_base::nrow (); i++)
			{
				for (j = 0; j < t_base::ncol (); j++)
					::printf ("%.0f \t", Element (i, j)) ;
				::printf ("\r\n") ;
			}
		}

		static IMat <T> diag (DWORD dwDim)
		{
			IMat<T>  ret (dwDim, dwDim) ;

			ret.Reset ((T) 0) ;

			DWORD i ;
			for (i = dwDim - 1; i != (DWORD) -1; i--)
				ret(i, i) = 1 ;
			return ret ;
		}

		IMat <T> setdiag () const
		{
			t_base::Reset () ;
			DWORD c = ::min (t_dims::ncol (), t_dims::nrow ()) - 1 ;
			for (;c != (DWORD) -1;c--)
				Element (c, c) = 1 ;
			return *this ;
		}

		IMat <T> setdiag2 () const
		{
			t_base::Reset () ;
			DWORD c = ::min (t_dims::ncol (), t_dims::nrow ()) - 1 ;
			DWORD cs = c ;

			for (;c != (DWORD) -1;c--)
				Element (cs - c, c) = 1 ;
			return *this ;
		}

		IMatByRow <T> &byrow () const { return * (IMatByRow <T> *) this ; }

		operator const IMatConst <T> & () const { return * (IMatConst <T> *) this ; }

		IMat <T> & Reshape (const DWORD r, const DWORD c)
		{
			if (r * c > t_base::m_pDataRef->GetSize ())
				return Recreate (r, c) ;
			IDim::ReCreate (r, 1, t_dims::m_apDim[0]) ;
			IDim::ReCreate (c, r, t_dims::m_apDim[1]) ;

			return *this ;
		}

		IMat <T> &Recreate (DWORD r, DWORD c)
		{
			t_base::Recreate (IIter <2> (r, c)) ;
			return *this ;
		}


	protected:

		inline T &Element (DWORD r, DWORD c) const { return t_base::m_pDataRef->GetItem (t_base::GetIdx (r, c) ) ; }
	} ;

	typedef IMat<double> IMatD ;
//	typedef IMatEdit<double> IMatEditD ;
	typedef IMatConst<double> IMatConstD ;
	typedef CSimpleArray <IMat <double> > IMatArrayD  ;
	typedef IMat<DWORD> IMatDW ;
//	typedef IMatEdit<DWORD> IMatEditDW ;
	typedef IMatConst<DWORD> IMatConstDW ;


///////////////////////////////////////////////////////////////////////////////
//	IMatByRow
///////////////////////////////////////////////////////////////////////////////

	template <class T>
	class IMatByRow
	{
		public:
			IMPL_OPERATOR_MAT_BY_ROW (+, +=, FC::FC_plus,		IMat, T)
			IMPL_OPERATOR_MAT_BY_ROW (/, /=, FC::FC_divide,		IMat, T)
			IMPL_OPERATOR_MAT_BY_ROW (*, *=, FC::FC_multiply,	IMat, T)
			IMPL_OPERATOR_MAT_BY_ROW (^, ^=, FC::FC_pow,		IMat, T)
//			IMPL_OPERATOR_MAT_BY_ROW (-, -=, FC::FC_minus,		IMat, T) ;

			template <class U> inline IMat<T>  operator -   (const IVec<U> &vec)	const	{ return FC_ElOp<FC::FC_minus, T>::OpMV_row	(Item (), vec); }
			template <class U> inline IMat<T> &operator -= (const IVec<U> &vec)		 	{ return FC_ElOpAs<FC::FC_minus>	::OpMV_row	(Item (), vec); }


		protected:

		inline IMat <T> &Item () const { return *(IMat <T> *) this ; }
	};

///////////////////////////////////////////////////////////////////////////////
//	IMatL
///////////////////////////////////////////////////////////////////////////////

	template <class T>
	class IMatL : public IMat <T>
	{
		public:
		IMPL_OPERATOR_VEC(|, |= , FC::FC_or,	IMatL, T)
		IMPL_OPERATOR_VEC(&, &= , FC::FC_and,	IMatL, T)
		IMPL_OPERATOR_VEC_NOAS(||, FC::FC_OR,	IMatL, bool)
		IMPL_OPERATOR_VEC_NOAS(&&, FC::FC_AND,	IMatL, bool)

		inline IVecL<bool> operator ! () const { return FC_ElOp<FC::FC_NOT, bool>::OpV (*this); }
		inline IVecL<T> operator ~ () const { return FC_ElOp<FC::FC_not, bool>::OpV (*this); }
	} ;

///////////////////////////////////////////////////////////////////////////////
//	Global matrix functions
///////////////////////////////////////////////////////////////////////////////

	template <class T> IMat <T> t (IMat <T> mat) { return mat.t () ; }


	template <class TA, class TB, class T>
	BOOL matmultmat (const TA &a, const TB &b, const IMat<T> &res)
	{
//		ASSERT (!a.IsLocked () && !b.IsLocked ()) ;
		ASSERT (a.ncol () == b.nrow ()) ;

		ASSERT (res.nrow () == a.nrow ()) ;
		ASSERT (res.ncol () == b.ncol ()) ;

		if (a.ncol () != b.nrow ())
			return FALSE ;		//	make an assertion!!! - check assertions!

		DWORD dwBncolm1 = b.ncol () - 1 ;
		DWORD dwAncolm1 = a.ncol () - 1 ;

		DWORD i, j, h ;
		for (i = a.nrow () - 1; i != (DWORD) -1; i--)
			for (j = dwBncolm1; j != (DWORD) -1; j--)
			{
				T &cur = res (i, j) = 0 ;
				for (h = dwAncolm1; h != (DWORD) -1; h--)
					 cur += (T) (a (i, h) * b(h, j)) ;
			}

		return TRUE ;
	}

	template <class TA, class TB, class TTC>
	BOOL vecmultmat (const TA &a, const TB &b, const IMat<TTC> &res)
	{
//		ASSERT (!a.IsLocked () && !b.IsLocked ()) ;
		ASSERT (a.size () == b.nrow ()) ;

		if (a.size () != b.nrow ())
			return FALSE ;		//	make an assertion!!! - check assertions!

		res.Create (1, b.ncol ()) ;

//		IMatEdit<TTC> resedit = res.Edit () ;

		res.Reset () ;

		DWORD dwAncolm1 = a.size () - 1 ;

		DWORD j, h ;
		for (j = b.ncol () - 1 ; j != (DWORD) -1; j--)
			for (h = dwAncolm1; h != (DWORD) -1; h--)
				res (0, j) += (TTC) (a (h) * b(h, j)) ;

		return TRUE ;
	}


	template <class TA, class TB, class TTC>
	BOOL matmultvec (const TA &a, const TB &b, const IMat<TTC> &res)
	{
//		ASSERT (!a.IsLocked () && !b.IsLocked ()) ;
		ASSERT (a.ncol () == b.size ()) ;
		ASSERT (res.size () == a.nrow ()) ;

		if (a.ncol () != b.size ())
			return FALSE ;		//	make an assertion!!! - check assertions!

//		res.Create (a.nrow (), 1) ;

//		IMatEdit<TTC> resedit = res.Edit () ;

		res.Reset () ;

		DWORD dwAncolm1 = a.ncol () - 1 ;

		DWORD i, h ;
		for (i = a.nrow () - 1; i != (DWORD) -1; i--)
				for (h = dwAncolm1; h != (DWORD) -1; h--)
					res (i, 0) += (TTC) (a (i, h) * b(h)) ;

		return TRUE ;
	}

	template <class TA, class TB, class TTC>
	BOOL matmultvec (const TA &a, const TB &b, const IVec<TTC> &res)
	{
//		ASSERT (!a.IsLocked () && !b.IsLocked ()) ;
		ASSERT (a.ncol () == b.size ()) ;
		ASSERT (res.size () == a.nrow ()) ;

		if (a.ncol () != b.size ())
			return FALSE ;		//	make an assertion!!! - check assertions!

//		res.Create (a.nrow (), 1) ;

//		IMatEdit<TTC> resedit = res.Edit () ;

		res.Reset () ;

		DWORD dwAncolm1 = a.ncol () - 1 ;

		DWORD i, h ;
		for (i = a.nrow () - 1; i != (DWORD) -1; i--)
				for (h = dwAncolm1; h != (DWORD) -1; h--)
					res (i) += (TTC) (a (i, h) * b(h)) ;

		return TRUE ;
	}


