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


#include "IMat.h"



//////////////////////////////////////////////////////////////////////////////
//	ISub_
//////////////////////////////////////////////////////////////////////////////

	BOOL ISub_::CreateIdx (IVec <DWORD> &adIdx, IDim *pdim) const
	{
		if (m_istType & ist_rem)
			return ((ISubR *)this)->CreateIdx (adIdx, pdim) ;
		return ((ISub *)this)->CreateIdx (adIdx, pdim) ;
	}

	ISub_::~ISub_ ()
	{
		if (m_istType & ist_del)
			if (m_ismMode == ism_ivecB)
				delete m_pbVec ;
			else if (m_ismMode == ism_ivecN || m_ismMode == ism_ivecNL)
				delete m_pdwVec ;
	}

//////////////////////////////////////////////////////////////////////////////
//	ISub	Intelligent Subindexing
///////////////////////////////////////////////////////////////////////////////

//	ISub::ISub (DWORD c1, DWORD c2, DWORD c3, ...):	ISub_ ((ISub_::ISubType) (ist_sel | ist_del), *new IVec<DWORD> (find (&c1, ARR_STOP))) {}
//	ISub::ISub (bool c1, ...) :						ISub_ ((ISub_::ISubType) (ist_sel | ist_del), *new IVec<bool> (find (&c1, ARR_STOP_B))) {}


	BOOL ISub::CreateIdx (IVec <DWORD> &adIdx, IDim *pdim) const
	{
		switch (m_ismMode)
		{
		case ism_none :		return FALSE ;				//	nothing 2do - did not change / create anything..
		case ism_idx :		return CreateIdx			(pdim, adIdx, m_dwIdx1) ;
		case ism_range :	return CreateIdx			(pdim, adIdx, m_dwIdx1, m_dwIdx2) ;
		case ism_ivecB :	return CreateIdx_VecBool	(pdim, adIdx, *m_pbVec) ;
		case ism_ivecN :	return CreateIdx_VecIdx		(pdim, adIdx, *m_pdwVec) ;
		case ism_ivecNL :	return CreateIdx_VecBool	(pdim, adIdx, *m_pdwVec) ;
		case ism_ptrB :		return CreateIdx			(pdim, adIdx, m_pbIdx, m_dwCount) ;
		case ism_ptrN :		return CreateIdx			(pdim, adIdx, m_pdwIdx, m_dwCount) ;
		case ism_ptrNL :	return CreateIdx_L			(pdim, adIdx, m_pdwIdx, m_dwCount) ;
		}
		return FALSE ;
	}

	BOOL ISub::CreateIdx (IDim *pdim, IVec <DWORD> &adIdx, DWORD dwIdx) const
	{
		if (dwIdx >= pdim->GetSize ())
			return FALSE ;
		adIdx.Create (1) ;
//		IVecEdit <DWORD> adIdx_edit () ;
		adIdx(0) = dwIdx ;
		return TRUE ;
	}

	BOOL ISub::CreateIdx (IDim *pdim, IVec <DWORD> &adIdx, DWORD dwIdxL, DWORD dwIdxU) const
	{
		if (dwIdxU > pdim->GetSize ())
			dwIdxU = pdim->GetSize () ;
		if (dwIdxL > dwIdxU)
			dwIdxL = dwIdxU ;

		adIdx.Create (dwIdxU - dwIdxL) ;

		DWORD i, v = 0 ; 

		for (i = dwIdxL; i < dwIdxU; i++)
			adIdx (v++) = pdim->Get_NC (i) ;

		return TRUE ;
	}

	template <class T>	//	T = DWORD or bool
	BOOL ISub::CreateIdx_VecBool (IDim *pdim, IVec <DWORD> &adIdx, const IVec<T> &vec) const
	{
		DWORD i, dwSize = 0 ;
/*		DWORD s1 = pdim->GetSize () ;
		DWORD d2 = vec.size () ;

		DWORD dwReadBool1 = min (1, 2) - 1 ;
*/
		DWORD dwReadBool = min (pdim->GetSize (), vec.size ()) - 1 ;
		for (i = dwReadBool; i != (DWORD) -1; i--)
			if (vec (i))
				dwSize ++ ;

		adIdx.Create (dwSize) ;

		for (i = dwReadBool; i != (DWORD) -1; i--)
			if (vec (i))
				adIdx  (--dwSize) = pdim->Get_NC (i) ;

		return TRUE ;
	}

//	template <class T>	//	T = IVec<DWORD> or IVec<DWORD>
	BOOL ISub::CreateIdx_VecIdx (IDim *pdim, IVec <DWORD> &adIdx, const IVec<DWORD> &vec) const
	{
		DWORD i, dwSize = 0 ;
		DWORD dwDimSize = pdim->GetSize () ;
		for (i = vec.size () - 1; i != (DWORD) -1; i--)
			if (vec(i) < dwDimSize)
				dwSize++ ;

		adIdx.Create (dwSize) ;
//		IVecEdit <DWORD> adIdx_edit (adIdx.Edit ()) ;

		for (i = vec.size () - 1; i != (DWORD) -1; i--)
		{
			DWORD dwCurIdx = vec(i) ;
			if (dwCurIdx < dwDimSize)
				adIdx (--dwSize) = pdim->Get_NC (dwCurIdx) ;
		}

		return TRUE ;
	}

	BOOL ISub::CreateIdx (IDim *pdim, IVec <DWORD> &adIdx, bool *pbIdx, DWORD dwCount) const
	{
		return CreateIdx_VecBool	(pdim, adIdx, IVec<bool> (dwCount, pbIdx)) ;
	}

	BOOL ISub::CreateIdx (IDim *pdim, IVec <DWORD> &adIdx, DWORD *pdwIdx, DWORD dwCount) const
	{
		return CreateIdx_VecIdx	(pdim, adIdx, IVec<DWORD> (dwCount, pdwIdx)) ;
	}

	BOOL ISub::CreateIdx_L (IDim *pdim, IVec <DWORD> &adIdx, DWORD *pdwIdx, DWORD dwCount) const
	{
		return CreateIdx_VecBool	(pdim, adIdx, IVec<DWORD> (dwCount, pdwIdx)) ;
	}


//////////////////////////////////////////////////////////////////////////////
//	ISubR	Intelligent Subindexing - removing elements
///////////////////////////////////////////////////////////////////////////////

	BOOL ISubR::CreateIdx (IVec <DWORD> &adIdx, IDim *pdim) const
	{
		switch (m_ismMode)
		{
		case ism_none :		return FALSE ;				//	nothing 2do - did not change / create anything..
		case ism_idx :		return CreateIdx			(pdim, adIdx, m_dwIdx1) ;
		case ism_range :	return CreateIdx			(pdim, adIdx, m_dwIdx1, m_dwIdx2) ;
		case ism_ivecB :	return CreateIdx_VecBool	(pdim, adIdx, *m_pbVec) ;
		case ism_ivecN :	return CreateIdx_VecIdx		(pdim, adIdx, *m_pdwVec) ;
		case ism_ivecNL :	return CreateIdx_VecBool	(pdim, adIdx, *m_pdwVec) ;
		case ism_ptrB :		return CreateIdx			(pdim, adIdx, m_pbIdx, m_dwCount) ;
		case ism_ptrN :		return CreateIdx			(pdim, adIdx, m_pdwIdx, m_dwCount) ;
		case ism_ptrNL :	return CreateIdx_L			(pdim, adIdx, m_pdwIdx, m_dwCount) ;
			
		}
		return FALSE ;
	}

	BOOL ISubR::CreateIdx (IDim *pdim, IVec <DWORD> &adIdx, DWORD dwIdx) const
	{
		DWORD i =  pdim->GetSize () - 1; 

		if (dwIdx > i)
			return FALSE ;

		DWORD v = i - 1 ;

		adIdx.Create (i) ;
//		IVecEdit <DWORD> adIdx_edit (adIdx.Edit ()) ;

		for (; i != (DWORD) -1; i--)
			if (i != dwIdx)
				adIdx(v--) = pdim->Get_NC (i) ;

		return TRUE ;
	}


	BOOL ISubR::CreateIdx (IDim *pdim, IVec <DWORD> &adIdx, DWORD dwIdxL, DWORD dwIdxU) const
	{
		if (dwIdxU > pdim->GetSize ())
			dwIdxU = pdim->GetSize () ;

		if (dwIdxL > dwIdxU)
			dwIdxL = dwIdxU ;

		DWORD dwSize = dwIdxL + pdim->GetSize () - dwIdxU ;

		adIdx.Create (dwSize) ;
//		IVecEdit <DWORD> adIdx_edit (adIdx.Edit ()) ;

		DWORD i ; 

		for (i = pdim->GetSize () - 1; i >= dwIdxU ; i--)
			adIdx (--dwSize) = pdim->Get_NC (i) ;

		for (i = dwIdxL - 1; i != (DWORD) -1 ; i--)
			adIdx (--dwSize) = pdim->Get_NC (i) ;

		return TRUE ;
	}

	template <class T>	//	T = DWORD or bool //IVec<bool> or IVec<bool>
	BOOL ISubR::CreateIdx_VecBool (IDim *pdim, IVec <DWORD> &adIdx, const IVec<T> &vec) const
	{
		DWORD i, dwSize = 0 ;
		DWORD dwReadBool = min (pdim->GetSize (), vec.size ()) - 1 ;
		for (i = dwReadBool; i != (DWORD) -1; i--)
			if (!vec (i))
				dwSize ++ ;

		adIdx.Create (dwSize) ;
//		IVecEdit <DWORD> adIdx_edit (adIdx.Edit ()) ;

		for (i = dwReadBool; i != (DWORD) -1; i--)
			if (!vec (i))
				adIdx  (--dwSize) = pdim->Get_NC (i) ;

		return TRUE ;
	}


//	template <class T>	//	T = IVec<DWORD> or IVec<DWORD>
	BOOL ISubR::CreateIdx_VecIdx (IDim *pdim, IVec <DWORD> &adIdx, const IVec<DWORD> &vec) const
	{
		DWORD i, dwDimSize = pdim->GetSize () ;
		bool *pb = new bool [dwDimSize] ;
		memset (pb, 1, dwDimSize) ;

		DWORD dwDel = 0 ;
		for (i = vec.size () - 1; i != (DWORD) -1; i--)
		{
			DWORD dwCurIdx = vec(i) ;
			if (dwCurIdx < dwDimSize && pb[dwCurIdx])
			{
				dwDel ++ ;
				pb[dwCurIdx] = false ;
			}
		}

		dwDel = dwDimSize - dwDel ;
		adIdx.Create (dwDel) ;
//		IVecEdit <DWORD> adIdx_edit (adIdx.Edit ()) ;

		for (i = dwDimSize - 1; i != (DWORD) -1; i--)
			if (pb[i])
				adIdx(--dwDel) = pdim->Get_NC (i) ;

		delete [] pb ;

		return TRUE ;
	}

	BOOL ISubR::CreateIdx (IDim *pdim, IVec <DWORD> &adIdx, bool *pbIdx, DWORD dwCount) const
	{
		return CreateIdx_VecBool (pdim, adIdx, IVec<bool> (dwCount, pbIdx)) ;
	}

	BOOL ISubR::CreateIdx (IDim *pdim, IVec <DWORD> &adIdx, DWORD *pdwIdx, DWORD dwCount) const
	{
		return CreateIdx_VecIdx (pdim, adIdx, IVec<DWORD> (dwCount, pdwIdx)) ;
	}

	BOOL ISubR::CreateIdx_L (IDim *pdim, IVec <DWORD> &adIdx, DWORD *pdwIdx, DWORD dwCount) const
	{
		return CreateIdx_VecBool	(pdim, adIdx, IVec<DWORD> (dwCount, pdwIdx)) ;
	}


//////////////////////////////////////////////////////////////////////////////
//	IDim	Intelligent Dimension Management
///////////////////////////////////////////////////////////////////////////////

	IDim::IDim (const IDim &idim) : m_bFull (idim.m_bFull)
	{
		m_dwSize = idim.m_dwSize ;
		if (m_bFull)
			m_dwStepSize = idim.m_dwStepSize ;
		else		//	this is a dimension with operational vector - copy this vector!
		{
			m_pdwSmartIDC = (DWORD *) malloc (m_dwSize * sizeof (DWORD)) ;
			memcpy (m_pdwSmartIDC, idim.m_pdwSmartIDC, m_dwSize * sizeof (DWORD)) ;
		}
	}

	IDim::~IDim ()
	{
		if (!m_bFull)
			::free (m_pdwSmartIDC) ;
	}

	void IDim::free ()
	{
		if (!m_bFull)
		{
			::free (m_pdwSmartIDC) ;
			m_pdwSmartIDC = 0 ;
		}
	}

	IDim *IDim::ReCreate (DWORD dwSize, DWORD dwStepSize, IDim *&pref)
	{
		if (pref)
		{
			if (pref->SingleOwner ())
			{
				pref->free () ;
				pref->m_bFull = true ;
				pref->m_dwSize = dwSize ;
				pref->m_dwStepSize = dwStepSize ;
				return pref ;
			}
			return Create (dwSize, dwStepSize, pref) ;
		}
		return Create (dwSize, dwStepSize, pref) ;
	}

	IDim *IDim::Create (const ISub_ &sub, IDim *pref)
	{
		if (sub.IsNone ())
			return pref ;

		IVec <DWORD> adwIdx ;
		if (!sub.CreateIdx (adwIdx, pref))
			return pref ;

		DWORD dwSize = adwIdx.size () ;
		return new IDim (0, dwSize, adwIdx.Detach ()) ;
	}

	IDim *IDim::Create (const ISub_ &sub, DWORD dwSize, DWORD dwStepSize)
	{
		IDim *pTemp = new IDim (dwSize, dwStepSize) ;

		if (sub.IsNone ())
			return pTemp ;

		IVec <DWORD> adwIdx ;
		if (!sub.CreateIdx (adwIdx, pTemp))
			return pTemp ;

		delete pTemp ;

		dwSize = adwIdx.size () ;
		return new IDim (0, dwSize, adwIdx.Detach ()) ;
	}


