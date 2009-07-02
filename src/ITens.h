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

#include "IDim.h"
#include "ElOp.h"
#include <R.h>


#define DIM_SEL		((DWORD) -1)
	
	template <class T> class IMat ;
//	template <class T> class IMatEdit ;
	template <class T> class IVec ;
//	template <class T> class IVecEdit ;
	template <class T, class D> class ITens_ ;
	template <class T, class D> class ITensFlat ;
//	template <class T, class D> class ITensEdit_ ;
	template <class T, class D> class ITensFlatEdit ;

	template <class T, class D>	class ITensConst ;

///////////////////////////////////////////////////////////////////////////////
//	IIter
///////////////////////////////////////////////////////////////////////////////
/*
	template <class T1> inline va_list GetArgList (const T1  *const pStart, ...)
	{
		va_list list ;
		va_start (list, (*pStart)) ;
		return list ;
	}
*/
/*	template <class T1> inline va_list GetArgList (const T1 &start, ...)
	{
		va_list list ;
		va_start (list, (start)) ;
		return list ;
	}
*/
#define OPERATOR_ITER(OP, OPAS, FUNC)	\
		inline IIter<D> operator OP		(const DWORD v) { IIter<D> ret ;	return * (IIter <D> *) FC_ElOp  <FC::FUNC, DWORD>::OpPE (m_adwIdx, D, v, ret.m_adwIdx) ; } \
		inline IIter<D> operator OPAS	(const DWORD v) {					return * (IIter <D> *) FC_ElOpAs<FC::FUNC>::OpPE (m_adwIdx, D, v) ; } \
		inline IIter<D> operator OP		(const IIter<D> &iter)	{ IIter<D> ret ; return * (IIter <D> *) FC_ElOp  <FC::FUNC, DWORD>::OpPP (m_adwIdx, D, iter.m_adwIdx, D, ret.m_adwIdx, D) ; } \
		inline IIter<D> operator OPAS	(const IIter<D> &iter)	{				 return * (IIter <D> *) FC_ElOpAs<FC::FUNC>::OpPP (m_adwIdx, D, iter.m_adwIdx, D) ; }

	template <DWORD D>
	class IIter 
	{
		public:

			IIter () {}

			IIter (DWORD c, ...)
			{
				va_list argList ;
				va_start (argList, c) ;

				DWORD d ;
				m_adwIdx[0] = c ;
				for (d = 1; d < D; d++)
					m_adwIdx[d] = va_arg(argList, DWORD) ;
			}

			inline IIter <D> &operator = (DWORD c)
			{
				DWORD d ;
				for (d = 0; d < D; d++)
					m_adwIdx[d] = c ;
				return *this ;
			}

			inline DWORD &operator [] (DWORD c) { return m_adwIdx [c] ; }
			inline const DWORD &operator [] (DWORD c) const { return m_adwIdx [c] ; }

//			inline operator DWORD * () { return m_adwIdx ; } 

			DWORD size () { DWORD d, ret = 0; for (d = 0; d < D; d++) ret *= m_adwIdx[d]; return d ; }
/*
		inline IIter<D> operator +	(const DWORD v) { IIter<D> ret ;	return * (IIter <D> *) FC_ElOp  <FC::FC_plus, DWORD>::OpPE (m_adwIdx, D, v, ret.m_adwIdx) ; }
		inline IIter<D> operator +=	(const DWORD v) {							return * (IIter <D> *) FC_ElOpAs<FC::FC_plus>::OpPE (m_adwIdx, D, v) ; }
		inline IIter<D> operator +	(const IIter<D> &iter)	{ IIter<D> ret ; return * (IIter <D> *) FC_ElOp  <FC::FC_plus, DWORD>::OpPP (m_adwIdx, D, iter.m_adwdx, D, ret.m_adwIdx, D) ; }
		inline IIter<D> operator +=	(const IIter<D> &iter)	{					return * (IIter <D> *) FC_ElOpAs<FC::FC_plus>::OpPP (m_adwIdx, D, iter.m_adwIdx, D) ; }
*/

			OPERATOR_ITER (+, +=, FC_plus)
			OPERATOR_ITER (-, -=, FC_minus)
			OPERATOR_ITER (*, *=, FC_multiply)
			OPERATOR_ITER (/, /=, FC_divide)
			OPERATOR_ITER (%, %=, FC_mod)

			bool mm (const IIter<D> &iter_max)
			{
				DWORD d ;
				for (d = 0; d < D; d++)
				{
					if (--m_adwIdx[d] != (DWORD) -1)
						return true ;
					m_adwIdx[d] = iter_max[d] - 1 ;
				}
				return false ;
			}

			const DWORD *GetIdxPtr () const { return m_adwIdx ; }


			const bool operator == (const IIter <D> &it)
			{
				return !memcmp (it.m_adwIdx, m_adwIdx, D * sizeof (DWORD)) ;
			}

			BOOL Dim2Iter (IIter <D> &it)
			{
				DWORD d ;
				for (d = 0; d < D; d++)
					if (m_adwIdx[d])
						it[d] = m_adwIdx[d] - 1 ;
					else
						return FALSE ;
				it[0] ++ ;
				return TRUE ;
			}
			
		protected:
			DWORD m_adwIdx [D] ;
	} ;

///////////////////////////////////////////////////////////////////////////////
//	IDimArray
///////////////////////////////////////////////////////////////////////////////

	typedef IDim * t_pIDim ;

	enum createMode
	{
		byRef,
		byVal
	} ;

	template <class T, DWORD D>
	class IDimArray
	{
		friend class ITens_<T, IDimArray <T, D> > ;
//		friend class ITensEdit_<T, IDimArray <T, D> > ;
		friend class ITensFlat<T, IDimArray <T, D> > ;
		

	protected:
			class SecType_1 {} ;
			class SecType_2 {} ;

			typedef SecType_1 t_cast ;
//			typedef SecType_2 t_editcast ;
	public:
		typedef IIter<D> t_Iter ;

		
		inline IDim *operator [] (DWORD d) const { return  m_apDim [d] ;}
		inline IDim &operator () (DWORD d) { return *m_apDim [d] ;}

		IDimArray () : m_dwOffset (0)
		{
			DWORD d ;
//Rprintf ("creating dimarray of size %d (this = 0x%08p) \r\n", D, this) ;
			for (d = 0; d < D; d++)
			{
//Rprintf ("setting dim %d to zero .. ", d) ;				
				m_apDim[d] = NULL ;
//Rprintf ("set \r\n", D) ;
//Rprintf ("emptying dim %d .. ", d) ;				
				//Empty (m_apDim[d] = NULL) ;
				IDim::Empty (m_apDim[d] = NULL) ;
//Rprintf ("emptyed \r\n", D) ;
			}
//Rprintf ("created dimarray of size %d\r\n", D) ;
		}

		IDimArray (const IDimArray<T, D> &ida, const createMode cm = byRef) : m_dwOffset (ida.m_dwOffset)
		{
			DWORD d ;
			if (cm == byRef)
				for (d = 0; d < D; d++)
					ida[d]->Reference (m_apDim[d] = NULL) ;
			else
			{
				DWORD dwCurStepSize = 1 ;
				for (d = 0; d < D; d++)
				{
					if (ida[d]->Full ())
						ida[d]->Reference (m_apDim[d] = NULL) ;
					else
						IDim::Create (ida[d]->GetSize (), dwCurStepSize, m_apDim[d] = NULL) ;
					dwCurStepSize *= m_apDim[d]->GetSize () ;
					m_dwOffset = 0 ;
				}
			}
		}

		IDimArray (const t_Iter &it) : m_dwOffset (0)
		{
			DWORD d, dwCurStepSize = 1 ;
			for (d = 0; d < D; d++)
			{
				IDim::ReCreate (it[d], dwCurStepSize, m_apDim[d] = NULL) ;
				dwCurStepSize *= it[d] ;
			}
		}

		IDimArray (DWORD c, va_list list) : m_dwOffset (0)
		{
			DWORD d, dwCurStepSize = 1 ;
			for (d = 0; d < D; d++)
			{
				if (d)
					c = va_arg (list, DWORD) ;
				IDim::ReCreate (c, dwCurStepSize, m_apDim[d] = NULL) ;
				dwCurStepSize *= c ;
			}

			//va_end (list) ;
		}

		IDimArray (DWORD dwOffset, IDim **ppDim) : m_dwOffset (dwOffset)
		{
			DWORD d ;
			for (d = 0; d < D; d++)
				ppDim[d]->Reference (m_apDim[d] = NULL) ;
		}

		IDimArray (DWORD dwOffset, IDim *pDim, ...) : m_dwOffset (dwOffset)
		{
			va_list list ;
			va_start (list, pDim) ;
			ConstructByArgList (pDim, list) ;
		}

		IDimArray (DWORD dwOffset, IDim *pDim, va_list list) : m_dwOffset (dwOffset)
		{
			ConstructByArgList (pDim, list) ;
		}

		IDimArray (const DWORD &dwSourceDims, const DWORD &dwOffset, IDim * const *ppDims, const ISub_ *pSub, va_list list) : m_dwOffset (dwOffset)
		{		//	dimension drop
			DWORD	dwCurIdx, d, cd = 0 ;//	current dim

			for (d = 0; d < dwSourceDims && cd < D; d++)
			{
				if (d)	pSub = &va_arg (list, ISub) ;
				if (pSub->IsIdx (dwCurIdx))
					m_dwOffset += ppDims[d]->Get_NC (dwCurIdx) ;
				else
					IDim::Create (*pSub, ppDims[d])->Reference (m_apDim[cd++] = NULL) ;
			}

//			if (cd >= D)	warning - too much dimensions selected

			for ( ;cd < D; cd++)
				Empty (m_apDim[cd++] = NULL) ;
		}

		IDimArray (const DWORD dwSourceDims, const DWORD dwOffset, IDim * const *ppDims, DWORD c, va_list list) : m_dwOffset (dwOffset)
		{		//	dimension drop

			DWORD d, cd = 0 ;	//	current dim

			for (d = 0; d < dwSourceDims; d++)
			{
				if (d) c = va_arg (list, DWORD) ;
				if (c == DIM_SEL)
					ppDims[d]->Reference (m_apDim[cd++] = NULL) ;
				else
					m_dwOffset += ppDims[d]->Get_NC (c) ;
			}

//			if (cd >= D)	warning - too much dimensions selected

			for ( ;cd < D; cd++)
				IDim::Empty (m_apDim[cd++] = NULL) ;
		}

		IDimArray (const DWORD dwSourceDims, const DWORD dwOffset, IDim * const *ppDims, const DWORD *pdwC) : m_dwOffset (dwOffset)
		{		//	dimension drop

			DWORD d, cd = 0 ;	//	current dim

			for (d = 0; d < dwSourceDims; d++)
			{
				//if (d) c = va_arg (list, DWORD) ;
				if (pdwC[d] == DIM_SEL)
					ppDims[d]->Reference (m_apDim[cd++] = NULL) ;
				else
					m_dwOffset += ppDims[d]->Get_NC (pdwC[d]) ;
			}

//			if (cd >= D)	warning - too much dimensions selected

			for ( ;cd < D; cd++)
				IDim::Empty (m_apDim[cd++] = NULL) ;
		}

		IDimArray (const ISub_ &s, va_list list, const IDimArray <T, D> &ida) : m_dwOffset (ida.m_dwOffset)
		{		//	subsetting dimensions
			DWORD d ;
			IDim::Create (s, ida[0])->Reference (m_apDim[0]) ;
			for (d = 1; d < D; d++)
				IDim::Create (va_arg (list, IDim), ida[d])->Reference (m_apDim[d]) ;
		}

		template <class T1>
		DWORD ReCreate (const IDimArray <T1, D> &it)
		{
			DWORD d, dwCurStepSize = 1 ;
			for (d = 0; d < D; d++)
			{
				
				IDim::ReCreate (it[d].size (), dwCurStepSize, m_apDim[d]) ;
				dwCurStepSize *= it[d].size() ;
			}
			return dwCurStepSize ;
		}

		DWORD ReCreate (const t_Iter &it)
		{
			DWORD d, dwCurStepSize = 1 ;
			for (d = 0; d < D; d++)
			{
				IDim::ReCreate (it[d], dwCurStepSize, m_apDim[d]) ;
				dwCurStepSize *= it[d] ;
			}
			return dwCurStepSize ;
		}

		DWORD ReCreate (DWORD c, va_list list)
		{

			IDim::ReCreate (c, 1, m_apDim[0]) ;
			DWORD d, dwCurStepSize = c ;
			for (d = 1; d < D; d++)
			{
				c = va_arg (list, DWORD) ;
				IDim::ReCreate (c, dwCurStepSize, m_apDim[d]) ;
				dwCurStepSize *= c ;
			}
			return dwCurStepSize ;
		}

		void ReCreate (DWORD dwOffset, IDim *pDim, va_list list)
		{
			m_dwOffset = dwOffset ;
			DWORD d ;
			pDim->Reference (m_apDim[0]) ;
			for (d = 1; d < D; d++)
				(pDim = va_arg (list, IDim *)) ->Reference (m_apDim[d]) ;
		}

/*
		void Reshape (const t_Iter &it)
		{
			DWORD d, dwCurStepSize = 1 ;
			for (d = 0; d < D; d++)
			{
				IDim::ReCreate (it[d], dwCurStepSize, m_apDim[d]) ;
				dwCurStepSize *= it[d] ;
			}
		}
*/
		BOOL IsFlat () const
		{
			DWORD d ;
			DWORD dwLastSize = 1 ;
			for (d = 0; d < D; d++)
			{
				IDim &dim = *m_apDim[d] ;
				if (!dim.Full  ())
					return FALSE ;

				if (dwLastSize != dim.GetStepSize ())
					return FALSE ;
				dwLastSize = dim.GetStepSize () * dim.GetSize () ;
			}
			return TRUE ;
		}

		~IDimArray ()
		{
//::Rprintf ("deleting dimarray 0x%08p of size %d\r\n", this, D) ;
			DWORD d ;
			for (d = 0; d < D; d++)
				CutRef (m_apDim[d]) ;
		}

		DWORD It_mm (t_Iter &it)
		{
			DWORD d ;
			for (d = 0; d < D; d++)
			{
				if (--it[d] != (DWORD) -1)
					return true ;
				it[d] = m_apDim[d]->GetSize () - 1 ;
			}
			return false ;
		}

		DWORD It_pp (t_Iter &it)
		{
			DWORD d ;
			for (d = 0; d < D; d++)
			{
				if (++it[d] < m_apDim[d]->GetSize ())
					return true ;
				it[d] = 0 ;
			}
			return false ;
		}

		inline DWORD size (DWORD c) const { return m_apDim[c]->GetSize () ; }

		DWORD size () const
		{
			register DWORD d, dwRet = 0 ;
			for (d = 0; d < D; d++)
				dwRet *= m_apDim [d]->GetSize () ;
			return dwRet ;
		}

//		DWORD GetIdx (DWORD c, ...) const { return GetIdx (*(t_Iter *) &c) ; }

		DWORD GetIdx (DWORD c, ...) const { va_list list ; va_start (list, c) ; return GetIdx (c, list) ; }

		DWORD GetIdx (DWORD c, va_list list) const
		{
			register DWORD idx = m_dwOffset + c, d ;
			for (d = 1; d < D; d++)
				idx += m_apDim[d]->Get_NC (/*c = */va_arg (list, DWORD)) ;
			//va_end (list) ;
			return idx ;
		}

		DWORD GetIdx (const t_Iter &it) const
		{
			register DWORD idx = m_dwOffset, d ;
			for (d = 0; d < D; d++)
				idx += m_apDim[d]->Get_NC (it[d]) ;
			return idx ;
		}

		t_Iter dim () const
		{
			t_Iter ret ;
			return dim (ret) ;
		}

		t_Iter &dim (t_Iter &it) const
		{
			DWORD d ;
			for (d = 0; d < D; d++)
				it[d] = m_apDim[d]->GetSize () ;
			return it ;
		}

		IDimArray <T, D> &Reference (IDimArray <T, D> &ida) const
		{
			DWORD d ;
			for (d = 0; d < D; d++)
				m_apDim[d]->Reference (ida.m_apDim [d]) ;
			ida.m_dwOffset = m_dwOffset ;
			return ida ;
		}

		inline DWORD ndim ()  const { return D ; }

//		void Reshape (DWORD dw, ...) { return 		}


	protected:

		void ConstructByArgList (IDim *pDim, va_list list)
		{
			DWORD d = 0 ;
			pDim->Reference (m_apDim[0] = NULL) ;
			Rprintf ("\tloading dim %d: %08p, size: %d\r\n", d, m_apDim[d], m_apDim[d]->GetSize ()) ;
			for (d = 1; d < D; d++)
			{
				pDim = va_arg (list, t_pIDim) ;
				pDim->Reference (m_apDim[d] = NULL) ;
				Rprintf ("\tloading dim %d: %08p, size: %d\r\n", d, pDim, 0) ; //m_apDim[d]->GetSize ()) ;
			}
		}

		DWORD m_dwOffset ;
		IDim *m_apDim[D] ;
	} ;

///////////////////////////////////////////////////////////////////////////////
//	IDimsMat
///////////////////////////////////////////////////////////////////////////////

	template <class T>
	class IDimsMat : public IDimArray <T, (DWORD) 2> 
	{
		friend class IMat <T> ;
//		friend class IMatEdit <T> ;
		friend class ITens_<T, IDimsMat <T> > ;
//		friend class ITensEdit_<T, IDimsMat <T> > ;
		friend class ITensFlat<T, IDimsMat <T> > ;
	protected:
		typedef IDimArray <T, (DWORD) 2> t_base ;
		typedef typename t_base::t_Iter t_Iter ;

//		typedef IMat<T> t_cast ;
//		typedef IMatEdit<T> t_editcast ;

	public:
		IDimsMat () : t_base () {}
		IDimsMat (const IDimsMat <T> &ida, const createMode cm = byRef) : t_base (ida, cm) {}
		IDimsMat (const t_Iter &it) : t_base (it) {}
			IDimsMat (DWORD c, va_list list) : t_base (c, list) {}
			
		IDimsMat (DWORD dwOffset, IDim **ppDim) : t_base (dwOffset, ppDim) {}
			IDimsMat (DWORD dwOffset, IDim *pDim, va_list list) : t_base (dwOffset, pDim, list) {}
		IDimsMat (DWORD dwSourceDims, DWORD dwOffset, IDim * const *ppDims, const ISub_ *pSub, va_list list) : t_base (dwSourceDims, dwOffset, ppDims, pSub, list) {}
		IDimsMat (DWORD dwSourceDims, DWORD dwOffset, IDim * const *ppDims, DWORD c, va_list list) : t_base (dwSourceDims, dwOffset, ppDims, c, list) {}
		IDimsMat (const DWORD dwSourceDims, const DWORD dwOffset, IDim * const *ppDims, const DWORD *pdwC) : t_base (dwSourceDims, dwOffset, ppDims, pdwC) {}
		IDimsMat (ISub_ *pSubs, const IDimsMat <T>  &ida) : t_base (pSubs, ida) {}

		inline DWORD nrow () const { return t_base::m_apDim[0]->GetSize () ; }
		inline DWORD ncol () const { return t_base::m_apDim[1]->GetSize () ; }

		inline DWORD size () const { return t_base::m_apDim[0]->GetSize () * t_base::m_apDim[1]->GetSize (); }
		inline DWORD size (DWORD d) const { return t_base::m_apDim[d]->GetSize () ; }

		typedef IMat <T> t_cast ;
//		typedef IMatEdit <T>  t_editcast ;

		DWORD It_mm (t_Iter &it)
		{
			if (--it[0] != (DWORD) -1)
				return 2 ;
			it[0] = t_base::m_apDim[0] - 1 ;

			if (--it[1] != (DWORD) -1)
				return 1 ;
			it[1] = t_base::m_apDim[1] - 1 ;
			return 0 ;
		}

		DWORD It_pp (t_Iter &it)
		{
			if (++it[0] < t_base::m_apDim[0])
				return 2 ;
			it[0] = 0 ;

			if (++it[1] < t_base::m_apDim[1])
				return 1 ;
			it[1] = 0 ;
			return 0 ;
		}

	protected:
		inline DWORD GetIdx (const t_Iter &it) const { return t_base::m_dwOffset + t_base::m_apDim [0]->Get_NC (it[0]) + t_base::m_apDim [1]->Get_NC (it[1]) ; } ;
		inline DWORD GetIdx (const DWORD  &r, const DWORD  &c) const { return t_base::m_dwOffset + t_base::m_apDim [0]->Get_NC (r) + t_base::m_apDim [1]->Get_NC (c) ; } ;
	} ;

	template <class T>
	class IDimsVec : public IDimArray <T, (DWORD) 1> 
	{
		friend class IVec <T> ;
//		friend class IVecEdit <T> ;
		friend class ITens_<T, IDimsVec <T> > ;
//		friend class ITensEdit_<T, IDimsVec <T> > ;
		friend class ITensFlat<T, IDimsVec <T> > ;
	protected:
		typedef IDimArray <T, (DWORD) 1> t_base ;
		typedef typename t_base::t_Iter t_Iter ;

//		typedef IVec<T> t_cast ;

	public:
		IDimsVec () : t_base () {}
		IDimsVec (const IDimsVec <T> &ida, const createMode cm = byRef) : t_base (ida, cm) {}
		IDimsVec (const t_Iter &it) : t_base (it) {}
			IDimsVec (DWORD c, va_list list) : t_base (c, list) {}
		IDimsVec (DWORD dwOffset, IDim **ppDim) : t_base (dwOffset, ppDim) {}
			IDimsVec (DWORD dwOffset, IDim *pDim, va_list list) : t_base (dwOffset, pDim, list) {}			
		IDimsVec (DWORD dwOffset, IDim *pDim) : t_base (dwOffset, pDim)	{}

		IDimsVec (DWORD dwSourceDims, DWORD dwOffset, IDim * const *ppDims, ISub_ *pSub, va_list list) : t_base (dwSourceDims, dwOffset, ppDims, pSub, list) {}
		IDimsVec (DWORD dwSourceDims, DWORD dwOffset, IDim * const *ppDims, DWORD c, va_list list) : t_base (dwSourceDims, dwOffset, ppDims, c, list) {}
		IDimsVec (const DWORD dwSourceDims, const DWORD dwOffset, IDim * const *ppDims, const DWORD *pdwC) : t_base (dwSourceDims, dwOffset, ppDims, pdwC) {}
		IDimsVec (ISub_ *pSubs, const IDimsVec <T> &ida) : t_base (pSubs, ida) {}

		inline DWORD size () const { return t_base::m_apDim[0]->GetSize () ; }
		inline DWORD size (DWORD d) const { return t_base::m_apDim[d]->GetSize () ; }

		typedef IVec <T> t_cast ;

		DWORD It_mm (t_Iter &it)
		{
			if (--it[0] != (DWORD) -1)
				return 1 ;
			it[0] = t_base::m_apDim[0] - 1 ;

			return 0 ;
		}

		DWORD It_pp (t_Iter &it)
		{
			if (++it[0] < t_base::m_apDim[0])
				return 1 ;
			it[0] = 0 ;

			return 0 ;
		}
	protected:
		inline DWORD GetIdx (const t_Iter &it)	const { return t_base::m_dwOffset + t_base::m_apDim [0]->Get_NC (it[0]) ; } ;
		inline DWORD GetIdx (const DWORD  &c)	const { return t_base::m_dwOffset + t_base::m_apDim [0]->Get_NC (c) ; } ;
	} ;

///////////////////////////////////////////////////////////////////////////////
//	ITens_
///////////////////////////////////////////////////////////////////////////////

#define ITens(TYPE, DIM)		ITens_		<TYPE, IDimArray<TYPE, (DWORD) DIM> >


class TClust ;
	template <class T, class D>
	class ITens_  : public D
	{
		friend class TClust ;
		friend class TClust1 ;
	protected:
		typedef D t_dims, t_base ;
		friend class ITensFlat <T, D> ;
		friend class ITensFlatEdit <T, D> ;

	public:
		typedef typename t_dims::t_Iter t_Iter ;
		typedef T m_datatype ;
		typedef ITensFlatEdit <T, D> t_flatedit	;
		typedef ITensFlat <T, D> t_flat ;
			// construction

//		ITens_ () { Empty (m_pDataRef = NULL) ; }
		ITens_ () { CDataRef<T>::Empty (m_pDataRef = NULL) ; }
		
		//ITens_ (DWORD c, ...) : t_dims (* (t_Iter*) &c) 	{ CDataRef<T>::Create (Dim(t_base::ndim () - 1)->GetFlatExtent (), m_pDataRef = NULL) ; }
		ITens_ (DWORD c, ...)// : t_dims (c, GetArgList (&c))
		{	
			va_list list ;
			va_start (list, c) ;
			t_dims::ReCreate (c, list) ;
			CDataRef<T>::Create (Dim(t_base::ndim () - 1)->GetFlatExtent (), m_pDataRef = NULL) ;
		}

		ITens_ (T *pData, SetDataFlag sdf, DWORD c, ...) //: t_base (* (t_Iter *) &c)
		{
			va_list list ;
			va_start (list, c) ;
			t_dims::ReCreate (c, list) ;
			CDataRef<T>::Create (Dim(t_dims::ndim () - 1)->GetFlatExtent (), m_pDataRef = NULL, pData, sdf, this) ;
		}
	//ITens_ (T *pData, SetDataFlag sdf, DWORD c, ...) : t_base (c, GetArgList (&c))	{	CDataRef<T>::Create (Dim(t_dims::ndim () - 1)->GetFlatExtent (), m_pDataRef = NULL, pData, sdf, this) ;	}


		//ITens_ (DWORD dwOffset, CDataRef <T> *pDataRef, IDim *pDim, ...) : t_dims (dwOffset, &pDim) //{ pDataRef->Reference (m_pDataRef = NULL) ; }
	ITens_ (DWORD dwOffset, CDataRef <T> *pDataRef, IDim *pDim, ...) //: t_dims (dwOffset, pDim, GetArgList (&pDim))
	{
			va_list list ;
			va_start (list, pDim) ;
			t_dims::ReCreate (dwOffset, pDim, list) ;

			pDataRef->Reference (m_pDataRef = NULL) ;
	}

		ITens_ (const t_Iter &it, T *pData = NULL, SetDataFlag sdf = sdf_modeReference) : D (it) { CDataRef<T>::Create (Dim(t_dims::ndim () - 1)->GetFlatExtent (), m_pDataRef = NULL, pData, sdf, this) ; }
		ITens_ (CDataRef <T> *pDataRef, const t_dims &ida) : t_dims (ida) { pDataRef->Reference (m_pDataRef = NULL) ; }

		ITens_ (const ITens_ <T, D> &tens, const createMode cm = byRef) : t_dims (tens, cm)
		{
			if (cm == byRef)
				tens.m_pDataRef->Reference (m_pDataRef = NULL) ;
			else
			{
				CDataRef <T>::Create (t_base::size (), m_pDataRef = NULL) ;
				FC_As::OpTT (*this, tens) ;
			}
		}

		~ITens_ ()
		{
//::Rprintf ("deleting ITens 0x%08p with dataref 0x%08p (dims:", this, m_pDataRef) ;

DWORD d ;
for (d = 0; d < t_dims::ndim (); d++)
//::Rprintf ("0x%08p-%d, ", t_dims::m_apDim[d], t_dims::m_apDim[d]->GetSize ()) ;

//::Rprintf (") ..\r\n") ;
			CutRef (m_pDataRef) ;
		}

		inline void Create (DWORD c, ...)
		{
			va_list list ;
			va_start (list, c) ;
			m_pDataRef->ReCreate (t_dims::ReCreate (c, list)) ;
		}

		void Create (const t_Iter &it, T *pData = NULL, SetDataFlag sdf = sdf_modeReference) 
		{ 
			m_pDataRef->ReCreate (t_dims::ReCreate (it), m_pDataRef, pData, sdf, this) ;
		}

		inline ITens_<T, D> operator () (const ISub_ &s, ...) { va_list list ; va_start (list, s) ; return ITens_<T, D> (m_pDataRef,	t_dims(s, list, *this)) ; }
		//inline ITens_<T, D> operator () (const ISub_ &s, va_list list) { return ITens_<T, D> (m_pDataRef,	t_dims(&s, list, *this)) ; }

		inline IVec <T> subIVec (const DWORD c, ...)			const { va_list list ; va_start (list, c) ; return subIVec (c, list) ; }
		inline IVec <T> subIVec (const DWORD c, va_list list)	const { return IVec<T> (m_pDataRef,	IDimsVec<T> (t_base::ndim (), t_dims::m_dwOffset, t_base::m_apDim, c, list)); }
		inline IVec <T> subIVec (const ISub_ &s, ...)			const {	va_list list ; va_start (list, s) ;return subIVec (s, list) ; }
		inline IVec <T> subIVec (const ISub_ &s, va_list list)	const {	return IVec<T> (m_pDataRef,	IDimsVec<T> (t_base::ndim (), t_dims::m_dwOffset, t_base::m_apDim, &s, list)); }
		inline IVec <T> subIVec (const t_Iter &it)				const { return IVec<T> (m_pDataRef,	IDimsVec<T> (t_base::ndim (), t_dims::m_dwOffset, t_base::m_apDim, it.GetIdxPtr ())); }

		inline IMat <T> subIMat (const DWORD c, ...)			const { va_list list ; va_start (list, c) ; return subIMat (c, list) ; }
		inline IMat <T> subIMat (const DWORD c, va_list list)	const { return IMat<T> (m_pDataRef,	IDimsMat<T> (t_base::ndim (), t_dims::m_dwOffset, t_base::m_apDim, c, list)); }
		inline IMat <T> subIMat (const ISub_ &s, ...)			const { va_list list ; va_start (list, s) ;return subIMat (s, list) ; }
		inline IMat <T> subIMat (const ISub_ &s, va_list list)	const { return IMat<T> (m_pDataRef,	IDimsMat<T> (t_base::ndim (), t_dims::m_dwOffset, t_base::m_apDim, &s, list)); }
		inline IMat <T> subIMat (const t_Iter &it)				const { return IMat<T> (m_pDataRef,	IDimsMat<T> (t_base::ndim (), t_dims::m_dwOffset, t_base::m_apDim, it.GetIdxPtr ())); }

/*		template<DWORD E>
		class sub
		{
		public:
			inline ITens_ <T, E> subITens (const DWORD c, ...)	const { return ITens_<T, E> (m_pDataRef, IDimArray <E> (D, IDimArray <T, D>::m_dwOffset, &IDimArray <T, D>::m_apDim, &c)); }
			inline ITens_ <T, E> subITens (const t_Iter &it)	const { return ITens_<T, E> (m_pDataRef, IDimArray <E> (D, IDimArray <T, D>::m_dwOffset, &IDimArray <T, D>::m_apDim, it)); }
			inline ITens_ <T, E> subITens (const ISub_ s, ...)	const { return ITens_<T, E> (m_pDataRef, IDimArray <E> (D, IDimArray <T, D>::m_dwOffset, &IDimArray <T, D>::m_apDim, &s)); }
		} ;
*/
		inline T &operator () (const DWORD c, ...)	const { va_list list ; va_start (list, c) ; return Element (c, list) ; }
		inline T &operator () (const t_Iter &it)	const { return Element (it) ;	}

		inline BOOL IsFlat () const { return t_dims::IsFlat () ; }

		inline ITensFlat<T, D> Flat () const	{ return ITensFlat<T, D> (*this) ; }
		inline ITensFlatEdit <T, D> FlatEdit ()	const { return ITensFlatEdit <T, D> (*this) ; }

		ITens_ <T, D> &operator = (const ITens_ <T, D> &tens)
		{
			tens.Reference (*this) ;
			tens.m_pDataRef->Reference (m_pDataRef) ;
			return *this ;
		}

		inline ITens_<T, D> CreateCopy () { return ITens_<T, D> (*this, byVal) ; }

		operator typename t_base::t_cast ()	 { return *(typename t_base::t_cast *) this ; }	//	cast to IMat or IVec

		inline void ToFlat (T *p) const { FC_As::OpRT (p, *this) ; }

		operator const ITensConst <T, D> & () const { return * (ITensConst <T, D> *) this ; }

		BOOL GetIMatArray (DWORD dwDimC, DWORD dwDimR, CSimpleArray <IMat <T> > & arr)
		{
			if (dwDimC == dwDimR || dwDimC >= t_base::ndim () || dwDimR >= t_base::ndim ())
				return FALSE ;

			DWORD d, dwSize = 1 ;
			t_Iter it, dim ;
			t_base::dim (dim) ;
			dim [dwDimR] = dim [dwDimC] = 1 ;

			it = dim - 1 ;

			for (d = t_base::ndim () - 1; d != (DWORD) -1 ; d--)
				if (d != dwDimR && d != dwDimC)
					dwSize *= t_base::m_apDim[d]->GetSize () ;

			arr.Create (dwSize) ;

			do
			{
				it [dwDimR] = it [dwDimC] = DIM_SEL ;
				arr[--dwSize] = subIMat (it) ;
				it [dwDimR] = it [dwDimC] = 0 ;
			} while (it.mm (dim)) ;

			return TRUE ;
		}

		BOOL GetIVecArray (DWORD dwDim, CSimpleArray <IVec <T> > & arr)
		{
			if (dwDim >= t_base::ndim ())
				return FALSE ;

			DWORD d, dwSize = 1 ;
			t_Iter it, dim ;
			t_base::dim (dim) ;
			dim [dwDim] = 1 ;

			it = dim - 1 ;

			for (d = t_base::ndim () - 1; d != (DWORD) -1 ; d--)
				if (d != dwDim)
					dwSize *= t_base::m_apDim[d]->GetSize () ;

			arr.Create (dwSize) ;

			do
			{
				it [dwDim] = DIM_SEL ;
				arr[--dwSize] = subIVec (it) ;
				it [dwDim] = 0 ;
			} while (it.mm (dim)) ;

			return TRUE ;
		}

		template <class U>
		const ITens_ <T, D> &operator << (const ITens_ <U, D> &tens) const
		{
			return FC_As::OpTT (*this, tens) ;
		}

		const ITens_ <T, D> &operator <<= (const ITens_ <T, D> &tens)
		{
			CDataRef <T>::Create (tens.size (), m_pDataRef) ;
			tens.Reference (*this) ;
			return FC_As::OpTT (*this, tens) ;
		}

		T max () const
		{
			if (!t_dims::size ())
				return 0 ;

			t_Iter	itDim = t_base::dim (),
					it = itDim - 1 ;

			T tMax = Element (it) ;

			while (it.mm (itDim))
			{
				T &dCur = Element (it) ;
				if (dCur > tMax) tMax = dCur ;
			}
			return tMax ;
		}

		BOOL AllEqual () const
		{
			if (!t_dims::size ())
				return FALSE ;

			t_Iter	itDim = t_base::dim (),
					it = itDim - 1 ;

			T &tFirst = Element (it) ;

			while (it.mm (itDim))
				if (tFirst != Element (it))
					return FALSE ;
			return TRUE ;
		}

		T sum () const
		{
			if (!t_dims::size ())
				return 0 ;

			t_Iter	itDim = t_base::dim (),
					it = itDim - 1 ;

			T tSum = Element (it) ;

			while (it.mm (itDim))
				tSum += Element (it) ;
			return tSum ;
		}

		T min () const
		{
			if (!t_dims::size ())
				return 0 ;

			t_Iter	itDim = t_base::dim (),
					it = itDim - 1 ;

			T tMin = Element (it) ;

			while (it.mm (itDim))
			{
				T &dCur = Element (it) ;
				if (dCur < tMin) tMin = dCur ;
			}
			return tMin ;
		}

		void MinMax (T &min, T&max) const
		{
			if (!t_dims::size ())
				return ;

			t_Iter	itDim = t_base::dim (),
					it = itDim - 1 ;

			min = Element (it) ;
			max = Element (it) ;

			while (it.mm (itDim))
			{
				T &dCur = Element (it) ;
				if (dCur > max)	max = dCur ;
				if (dCur < min)	min = dCur ;
			}
			
		}

		void Limit (T min, T max) const
		{
			if (!t_dims::size ())
				return ;
			if (min == max)
			{
				Reset (min) ;
				return ;
			}

			t_Iter	itDim = t_base::dim (),
					it = itDim - 1 ;

			if (min > max)
				swap (&min, &max) ;

			do
			{
				T &dCur = Element (it) ;
				if (dCur > max)
					dCur = max ;
				else if (dCur < min)
					dCur = min ;
			}
			while (it.mm (itDim)) ;
		}

		void Limit_L (T min) const
		{
			if (!t_dims::size ())
				return ;

			t_Iter	itDim = t_base::dim (),
					it = itDim - 1 ;

			do
			{
				T &dCur = Element (it) ;
				if (dCur < min) dCur = min ;
			}
			while (it.mm (itDim)) ;
		}

		void Limit_U (T max) const
		{
			if (!t_dims::size ())
				return ;

			t_Iter	itDim = t_base::dim (),
					it = itDim - 1 ;

			do
			{
				T &dCur = Element (it) ;
				if (dCur > max) dCur = max ;
			}
			while (it.mm (itDim)) ;
		}

		template <class U> inline const ITens_ <T, D> &operator += (const U &v) const { return FC_ElOpAs<FC::FC_plus>::OpTE (*this, v); }
		template <class U> inline const ITens_ <T, D> &operator -= (const U &v) const { return FC_ElOpAs<FC::FC_minus>::OpTE (*this, v); }
		template <class U> inline const ITens_ <T, D> &operator *= (const U &v) const { return FC_ElOpAs<FC::FC_multiply>::OpTE (*this, v); }
		template <class U> inline const ITens_ <T, D> &operator /= (const U &v) const { return FC_ElOpAs<FC::FC_divide>::OpTE (*this, v); }

		const ITens_<T, D> &operator = (const T &v) const { FC_As::OpTE (*this, v) ; }
		void Reset (const T  &val = 0) const { FC_As::OpTE (*this, val) ; }

//		template <class T>
		DWORD CountGreater (T val) const
		{
			if (!t_dims::size ())
				return 0 ;

			DWORD dwCount = 0 ;
			t_Iter	itDim = t_base::dim (),
					it = itDim - 1 ;

			do
			{
				if( Element (it) > val)
					dwCount ++ ;
			}
			while (it.mm (itDim)) ;
			return dwCount ;
		}

/*		DWORD CountGetGreater (const T &val, ITens_<bool, IDimArray<bool, D> &bRet) const
		{
			if (!t_dims::size ())
				return 0 ;

			bRet.Reshape (*this) ;

			DWORD dwCount = 0 ;
			t_Iter	itDim = t_base::dim (),
					it = itDim - 1 ;

			do
			{
				if( bRet(itDim) = Element (it) > val)
					dwCount ++ ;
			}
			while (it.mm (itDim)) ;
			return dwCount ;
		}
*/

/*
		template <class FUNC>
		DWORD count (const FUNC &func)
		{
			if (!t_dims::size ())
				return 0 ;

			DWORD dwCount = 0 ;
			t_Iter	itDim = t_base::dim (),
					it = itDim - 1 ;
			do
			{
				if (func (Element (it)))
					dwCount ++ ;
			}
			while (it.mm (itDim)) ;
			return dwCount ;
		}*/

		DWORD CountLess (T val) const
		{
			if (!t_dims::size ())
				return 0 ;

			DWORD dwCount = 0 ;
			t_Iter	itDim = t_base::dim (),
					it = itDim - 1 ;

			do
			{
				if( Element (it) < val)
					dwCount ++ ;
			}
			while (it.mm (itDim)) ;
			return dwCount ;
		}

/*		DWORD CountOuter (T val1, T Val2) const
		{
			if (!t_dims::size ())
				return 0 ;

			DWORD dwCount = 0 ;
			t_Iter	itDim = t_base::dim (),
					it = itDim - 1 ;

			do
			{
				T &dCur = Element (it) ;
				if( dCur < val1 || dCur > val2)
					dwCount ++ ;
			}
			while (it.mm (itDim)) ;
			return dwCount ;
		}
*/
		bool Equal (const ITens_ <T, D> &t) const
		{
			t_Iter	itDim = t_base::dim (),
					it = itDim - 1 ;

			if (t.dim () == itDim)
			{
				if (!itDim.size ())
					return true ;
				do
				{
					if (Element (it) != t(it))
						return false ;
				} while (it.mm (itDim)) ;
			}
			else
			{
				if (!itDim.size ())
					return false ;

				t_Iter	itDim2 = t.dim (),
						it2 = it % itDim2 ;

				do
				{
					if (Element (it) != t(it2))
						return false ;
					it2.mm (itDim2) ;
				} while(it.mm (itDim)) ;
			}
			return true ;
		}

		ITens_<T, D> & Reshape (const t_Iter &it)
		{
			if (it.size () > m_pDataRef->GetSize ())
				return Recreate (it) ;
			t_base::Reshape (it) ;
			return *this ;
		}

		template <class T1>
		ITens_<T, D> & Reshape (const ITens_<T1, D> &it)
		{
			if (it.size () > m_pDataRef->GetSize ())
				return Recreate (it) ;
			t_base::Reshape (it) ;
			return *this ;
		}

		ITens_<T, D> & Recreate (const t_Iter &it)
		{
			CDataRef<T>::Create (t_dims::ReCreate (it), m_pDataRef) ;
			return *this ;
		}

		template <class T1>
		ITens_<T, D> & Recreate (const ITens_<T1, D> &it)
		{
			CDataRef<T>::Create (t_dims::ReCreate (it), m_pDataRef) ;
			return *this ;
		}

	protected:

		inline T *Data () const { return m_pDataRef->GetData () ; }

		inline DWORD GetRefCount () const { return m_pDataRef->GetRefCount () ; }

		inline void MakeLockable ()
		{
			if (m_pDataRef->GetRefCount () > 1 && m_pDataRef->GetOwner () == this)
				CDataRef<T>::Create (m_pDataRef->GetSize (), m_pDataRef, m_pDataRef->CreateCopyDetach (), sdf_modeReference) ;
		}

//		inline T &Element (const t_Iter &it) { return m_pDataRef->GetData ()[GetIdx (it)] ; }

		inline T &Element (DWORD c, ...) const { va_list list ; va_start (list, c) ; return Element (c, list) ; }
		inline T &Element (DWORD c, va_list list) const { return m_pDataRef->GetItem (t_base::GetIdx (c, list)) ; }
		inline T &Element (const t_Iter &it) const { return m_pDataRef->GetItem (GetIdx (it)) ; }

		inline IDim *&Dim (const DWORD d) { return t_dims::m_apDim [d] ; }
		inline IDim *Dim (const DWORD d) const { return t_dims::m_apDim [d] ; }

		CDataRef<T> *m_pDataRef ;
	} ;

///////////////////////////////////////////////////////////////////////////////
//	ITensFlat
///////////////////////////////////////////////////////////////////////////////

	template <class T, class D>
	class ITensFlat
	{
		public:
			ITensFlat (const ITens_ <T, D> &tens)
			{
				if (tens.IsFlat ())
				{
					tens.m_pDataRef->Reference (m_pDataRef = NULL) ;
					m_dwOffset = tens.m_dwOffset ;
				}
				else
				{
					m_dwOffset = 0 ;
					CDataRef <T>::Create (tens.size (), m_pDataRef = NULL) ;

					FC_As::OpRT (m_pDataRef->GetData (), tens) ;
				}
			}

			~ITensFlat () { CutRef (m_pDataRef) ; }

			inline operator const T *	() { return m_pDataRef->GetData () + m_dwOffset ; }
//			inline T * HackConst		() { return m_pDataRef->GetData () + m_dwOffset ; }

		protected:
			CDataRef<T> *m_pDataRef ;
			DWORD m_dwOffset ;
	} ;

///////////////////////////////////////////////////////////////////////////////
//	ITensFlatEdit
///////////////////////////////////////////////////////////////////////////////

	template <class T, class D>
	class ITensFlatEdit
	{
		public:

			ITensFlatEdit (const ITens_ <T, D> &tens) : m_pTens (NULL)
			{
				m_pdwRefCount = new DWORD ;
				*m_pdwRefCount = 1 ;

				if (tens.IsFlat ())
				{
					tens.m_pDataRef->SoftLock (m_pDataRef = NULL) ;
					m_dwOffset = tens.m_dwOffset ;
				}
				else
				{
					m_pTens = &tens ;
					CDataRef <T>::Create (tens.size (), m_pDataRef = NULL) ;
					FC_As::OpRT (m_pDataRef->GetData (), tens) ;
					m_dwOffset = 0 ;
				}
			}

			ITensFlatEdit (const ITensFlatEdit <T, D> &tfe) : m_pDataRef (tfe.m_pDataRef), m_pTens (tfe.m_pTens), m_dwOffset (tfe.m_dwOffset), m_pdwRefCount (tfe.m_pdwRefCount)
			{
				*m_pdwRefCount++ ;
			}

			~ITensFlatEdit ()
			{
				if (--*m_pdwRefCount)
					return ;

				if (m_pTens)
				{
					FC_As::OpTR (*m_pTens, m_pDataRef->GetData ()) ;
					CutRef (m_pDataRef) ;
				}
				else
					Unlock (m_pDataRef) ;
				delete m_pdwRefCount ;
			}

			inline operator T *	() const { return m_pDataRef->GetData () + m_dwOffset ; }

		protected:
			CDataRef<T> *m_pDataRef ;
			const ITens_ <T, D> *m_pTens ;
			DWORD m_dwOffset ;
			DWORD *m_pdwRefCount ;
	} ;

///////////////////////////////////////////////////////////////////////////////
//	ITensConst
///////////////////////////////////////////////////////////////////////////////

	template <class T, class D>
	class ITensConst : public D
	{
	public:
		typedef typename D::t_Iter t_Iter ;

		inline T & operator ()	(DWORD c, ...)		const { va_list list ; va_start (list, c) ; return Element (c, list) ; }	//{ return Element (GetIdx (* (t_Iter *)&c)) ; }
		inline T & operator ()	(const t_Iter &it)	const { return Element (it) ; }				//{ return Element (GetIdx (it)) ; }

		inline T & Element		(DWORD c, ...) const	{ return Element (* (t_Iter *)&c) ; }
		inline T & Element		(const t_Iter &it) const { return m_pDataRef->GetData () [D::GetIdx (it)] ; }

		inline ITensFlat<T, D> Flat () const { return ITensFlat<T, D> (*(const ITens_ <T, D> *)this) ; }

	protected:
		CDataRef <T> *m_pDataRef ;
	} ;


