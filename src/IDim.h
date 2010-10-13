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


#include "IRef.h"

#define ARR_STOP	((DWORD) -1)
#define STOP		((DWORD) -1)
#define ARR_STOP_B	( *(bool *)"\xff")

	class IDim ;

	template <class T> class IMat ;
	template <class T> class IVec ;

///////////////////////////////////////////////////////////////////////////////
//	ISub	implements a simple subset class, which is then interpreted
///////////////////////////////////////////////////////////////////////////////

	class ISub_
	{
		friend class IDim ;
	protected:
		enum ISubMode { ism_none, ism_idx, ism_range, ism_ivecB, ism_ivecN, ism_ivecNL, ism_ptrB, ism_ptrN, ism_ptrNL } ;
		enum ISubType { ist_sel = 0x0, ist_rem = 0x1, ist_del = 0x2 } ;

	public:
		ISub_ (ISubType ist) :												m_ismMode (ism_none),	m_istType (ist) {}
		ISub_ (ISubType ist, DWORD dwIdx) :									m_ismMode (ism_idx),	m_istType (ist), m_dwIdx1 (dwIdx) {}
		ISub_ (ISubType ist, DWORD dwIdxL, DWORD dwIdxU) :					m_ismMode (ism_range),	m_istType (ist), m_dwIdx1 (dwIdxL), m_dwIdx2 (dwIdxU) {}

		ISub_ (ISubType ist, const IVec <bool>	& vec) :					m_ismMode (ism_ivecB),	m_istType (ist), m_pbVec (& vec)	{}
		ISub_ (ISubType ist, const IVec <DWORD>	& vec, BOOL bLogical) :		m_ismMode (bLogical ? ism_ivecNL : ism_ivecN),	m_istType (ist), m_pdwVec (& vec)	{}

		ISub_ (ISubType ist, bool *pbIdx, DWORD dwCount) :					m_ismMode (ism_ptrB),	m_istType (ist), m_pbIdx (pbIdx), m_dwCount (dwCount) {}
		ISub_ (ISubType ist, DWORD *pdwIdx, DWORD dwCount, BOOL bLogical) :	m_ismMode (bLogical ? ism_ptrNL : ism_ptrN), m_istType (ist), m_pdwIdx (pdwIdx), m_dwCount (dwCount) {}

		~ISub_ () ;

		inline bool IsNone () const { return m_ismMode == ism_none ; }
		inline bool IsIdx () const { return m_ismMode == ism_idx ; }
		inline bool IsIdx (DWORD &dwIdx) const
		{
			if (m_ismMode != ism_idx)
				return false ;
			dwIdx = m_dwIdx1;
			return true ;
		}

		BOOL CreateIdx (IVec <DWORD> &adIdx, IDim *pdim) const ;

	protected:

		ISubMode m_ismMode ;
		ISubType m_istType ;

		union { public:
			const IVec <bool>		*m_pbVec ;
			const IVec <DWORD>	*m_pdwVec ;
			bool			*m_pbIdx ;
			DWORD			*m_pdwIdx ;
			DWORD			m_dwIdx1 ;
		} ;

		union { public:
			DWORD	m_dwIdx2 ;
			DWORD	m_dwCount ;
		} ;
	} ;


	class ISub : public ISub_
	{
		friend class TClust ;
		friend class ISub_ ;
	public:
		ISub () :									ISub_ (ist_sel) {}
		ISub (DWORD dwIdx) :						ISub_ (ist_sel, dwIdx) {}
		ISub (DWORD dwIdxL, DWORD dwIdxU) :			ISub_ (ist_sel, dwIdxL, dwIdxU) {}

		ISub (const IVec <bool>	& vec) :			ISub_ (ist_sel, vec) {}
		ISub (const IVec <DWORD>& vec, BOOL bLogical = FALSE) :	ISub_ (ist_sel, vec, bLogical) {}

		ISub (bool *pbIdx, DWORD dwCount) :			ISub_ (ist_sel, pbIdx, dwCount) {}
		ISub (DWORD *pdwIdx, DWORD dwCount, BOOL bLogical = FALSE) :	ISub_ (ist_sel, pdwIdx, dwCount, bLogical) {}
//		ISub (DWORD c1, DWORD c2, DWORD c3, ...) ;
//		ISub (bool c1, ...) ;

	protected:
		BOOL CreateIdx (IVec <DWORD> &adIdx, IDim *pdim) const ;

		template <class T>
		BOOL CreateIdx_VecBool	(IDim *pdim, IVec <DWORD> &adIdx, const IVec <T> &vec) const ;
		BOOL CreateIdx_VecIdx	(IDim *pdim, IVec <DWORD> &adIdx, const IVec <DWORD> & vec) const ;

		BOOL CreateIdx (IDim *pdim, IVec <DWORD> &adIdx, DWORD dwIdx) const ;
		BOOL CreateIdx (IDim *pdim, IVec <DWORD> &adIdx, DWORD dwIdxL, DWORD dwIdxU) const ;

		BOOL CreateIdx (IDim *pdim, IVec <DWORD> &adIdx, bool *pbIdx, DWORD dwCount) const ;
		BOOL CreateIdx (IDim *pdim, IVec <DWORD> &adIdx, DWORD *pdwIdx, DWORD dwCount) const ;
		BOOL CreateIdx_L (IDim *pdim, IVec <DWORD> &adIdx, DWORD *pdwIdx, DWORD dwCount) const ;
	} ;

	class ISubR : public ISub_
	{
		friend class ISub_ ;
	public:
		ISubR () :								ISub_ (ist_rem) {}
		ISubR (DWORD dwIdx) :					ISub_ (ist_rem, dwIdx) {}
		ISubR (DWORD dwIdxL, DWORD dwIdxU) :	ISub_ (ist_rem, dwIdxL,	dwIdxU) {}

		ISubR (const IVec <bool>	& vec) :			ISub_ (ist_rem, vec) {}
		ISubR (const IVec <DWORD>	& vec, BOOL bLogical = FALSE) :			ISub_ (ist_rem, vec, bLogical) {}

		ISubR (bool *pbIdx, DWORD dwCount) :	ISub_ (ist_rem, pbIdx,	dwCount) {}
		ISubR (DWORD *pdwIdx, DWORD dwCount, BOOL bLogical = FALSE) :	ISub_ (ist_rem, pdwIdx,	dwCount, bLogical) {}

	protected:

		BOOL CreateIdx (IVec <DWORD> &adIdx, IDim *pdim) const ;

		BOOL CreateIdx (IDim *pdim, IVec <DWORD> &adIdx, DWORD dwIdx) const ;
		BOOL CreateIdx (IDim *pdim, IVec <DWORD> &adIdx, DWORD dwIdxL, DWORD dwIdxU) const ;

		inline BOOL CreateIdx (IDim *pdim, IVec <DWORD> &adIdx, bool *pbIdx, DWORD dwCount) const ;
		inline BOOL CreateIdx (IDim *pdim, IVec <DWORD> &adIdx, DWORD *pdwIdx, DWORD dwCount) const ;

		template <class T>
		BOOL CreateIdx_VecBool	(IDim *pdim, IVec <DWORD> &adIdx, const IVec <T> &vec) const ;
		BOOL CreateIdx_VecIdx	(IDim *pdim, IVec <DWORD> &adIdx, const IVec <DWORD> & vec) const ;
		BOOL CreateIdx_L (IDim *pdim, IVec <DWORD> &adIdx, DWORD *pdwIdx, DWORD dwCount) const ;
	} ;

///////////////////////////////////////////////////////////////////////////////
//	IDim	Intelligent Dimension Management
///////////////////////////////////////////////////////////////////////////////

//	class IDimEdit ;

	class TClust ;
	class IDim: public CRefBase <IDim>
	{
		friend class TClust ;
		friend class TClust1 ;
		friend class ISub ;
//		friend class IDimEdit ;
	public:

		DEC_EMPTY (IDim) ;

		IDim (DWORD dwSize = 0, DWORD dwStepSize = 1) : CRefBase <IDim> (1), m_dwSize (dwSize) , m_dwStepSize (dwStepSize), m_bFull (true) {}
				//	this is the public constructor, which is used when an object of this class is constructed "in the usual way" - since an object like that shall not be destroyed automatically, the refcount is set to 1

		IDim (const IDim &idim) ;
		static inline IDim *Create (DWORD dwSize, DWORD dwStepSize, IDim *&pIDim)
		{
			return (new IDim (0, dwSize, dwStepSize))->Reference (pIDim) ;
		}

		static IDim *Create (const ISub_ &sub, IDim *pref) ;
		static IDim *Create (const ISub_ &sub, DWORD dwSize, DWORD dwStepSize) ;

		static IDim *ReCreate (DWORD dwSize, DWORD dwStepSize, IDim *&pref) ;

		~IDim () ;

		void free () ;

		inline DWORD Get_NC (DWORD dwIDX) const
		{
			ASSERT (dwIDX < m_dwSize) ;
			switch (m_bFull)
			{
			case 0 : return m_pdwSmartIDC [dwIDX] ;
			case 1 : return dwIDX * m_dwStepSize ;
			}
			return  0 ;
		}

		inline DWORD GetSize ()		const { return m_dwSize ; }
		inline DWORD GetStepSize ()	const { return m_dwStepSize ; }
		inline BOOL	 Full ()		const { return m_bFull ; }

		inline DWORD GetFlatExtent () const { return m_dwSize * m_dwStepSize ; }

	protected:

		IDim (DWORD dwRefCount, DWORD dwSize, DWORD *pdwSmartIDC)	: CRefBase <IDim> (dwRefCount), m_dwSize (dwSize) , m_pdwSmartIDC (pdwSmartIDC), m_bFull (false) {}
		IDim (DWORD dwRefCount, DWORD dwSize, DWORD dwStepSize)		: CRefBase <IDim> (dwRefCount), m_dwSize (dwSize) , m_dwStepSize (dwStepSize), m_bFull (true) {}

		IDim *ReCreate_m (DWORD dwSize, DWORD dwStepSize, IDim *&pref) ;

		DWORD	m_dwSize ;
		union { public:
			DWORD	m_dwStepSize,	*m_pdwSmartIDC ;
		} ;
		bool	m_bFull ;
	} ;

	typedef CRefRef <IDim> IDimRef ;


