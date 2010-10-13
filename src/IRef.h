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

#include <memory.h>
//#include <stdlib.h>
#include <cstdarg>

#include <R.h>

typedef unsigned int DWORD ;
typedef int	BOOL ;
typedef unsigned char BYTE ;


#ifndef LPCTSTR
#define LPCTSTR const char *
#endif

#define protected public

#define TRUE	1
#define FALSE	0
#ifndef NULL
	#define NULL	0
#endif


#define BACK_CAST(CLASS, TO, FROM) ((TO *)(1 + ((BYTE *) CLASS) - ((BYTE *) (FROM *) (TO *) 1)))

#ifndef ASSERT

	#ifdef _DEBUG
		inline void  assertfunc (int b) { if (!b) b = * (int *) 0 ;}
		#define ASSERT(foo) //assertfunc (foo)
	#else
		#define ASSERT(FOO)
	#endif

#endif

///////////////////////////////////////////////////////////////////////////////
//	auxiliary functions / classes
///////////////////////////////////////////////////////////////////////////////

	template <class T>
	inline void swap (T *p1, T *p2)
	{
		BYTE abTemp [sizeof (T)] ;
		memcpy (abTemp, p1, sizeof (T)) ;
		memcpy (p1, p2, sizeof (T)) ;
		memcpy (p2, abTemp, sizeof (T)) ;
	}

	template <class T>
	DWORD find (T *p, const T &v)
	{
		DWORD c ;
		for (c = 0; *p != v; c++)
			p++ ;
		return c ;
	}

	template <class T>
	DWORD count (const T *p, DWORD dwSize, const T &v)
	{
		DWORD dwCount = 0 ;
		for ( ; dwSize; dwSize--)
			if (*p++ == v)
				dwCount++ ;
		return dwCount ;
	}

	template <class T> inline T sqr (const T &v) { return v * v ; }
		

	template <class T> inline const T &max (const T &a, const T &b) { return a > b ? a : b ; }
	template <class T> inline const T &min (const T &a, const T &b) { return a > b ? b : a ; }

//	template <class T> inline T &max (T &a, T &b) { return a > b ? a : b ; }
//	template <class T> inline T &min (T &a, T &b) { return a > b ? b : a ; }


	class FC
	{
	public:
		class FC_plus		{ public :	template <class TA, class TB, class TC> inline static void Calc (const TA &a,const TB &b, TC &c) { c = (TC) (a + b) ; }	} ;
		class FC_minus		{ public :	template <class TA, class TB, class TC> inline static void Calc (const TA &a,const TB &b, TC &c) { c = (TC) (a - b) ; }	} ;
		class FC_divide		{ public :	template <class TA, class TB, class TC> inline static void Calc (const TA &a,const TB &b, TC &c) { c = (TC) (a / b) ; }	} ;
		class FC_multiply	{ public :	template <class TA, class TB, class TC> inline static void Calc (const TA &a,const TB &b, TC &c) { c = (TC) (a * b) ; }	} ;
		class FC_pow		{ public :	template <class TA, class TB, class TC> inline static void Calc (const TA &a,const TB &b, TC &c) { c = (TC) pow (a, b) ; }} ;

		class FC_mod		{ public :	template <class TA, class TB, class TC> inline static void Calc (const TA &a,const TB &b, TC &c) { c = (TC) (a % b) ; }	} ;
		class FC_minus_u	{ public :	template <class T> inline static void Calc (T a, T &b) { b = -a ; }	} ;

		class FC_exp		{ public :	template <class TA, class TB> inline static void Calc (const TA &a, TB &b) { b = (TB) exp ((double) a) ; }	} ;
		class FC_log		{ public :	template <class TA, class TB> inline static void Calc (const TA &a, TB &b) { b = (TB) log ((double) a) ; }	} ;
		class FC_sqr		{ public :	template <class TA, class TB> inline static void Calc (const TA &a, TB &b) { b = (TB) (a * a) ; }	} ;

		class FC_or			{ public :	template <class TA, class TB, class TC> inline static void Calc (const TA &a,const TB &b, TC &c) { c = (TC) (a | b) ; }	} ;
		class FC_and		{ public :	template <class TA, class TB, class TC> inline static void Calc (const TA &a,const TB &b, TC &c) { c = (TC) (a & b) ; }	} ;
		class FC_OR			{ public :	template <class TA, class TB, class TC> inline static void Calc (const TA &a,const TB &b, TC &c) { c = (TC) (a || b) ; }	} ;
		class FC_AND		{ public :	template <class TA, class TB, class TC> inline static void Calc (const TA &a,const TB &b, TC &c) { c = (TC) (a && b) ; }	} ;
		class FC_xor		{ public :	template <class TA, class TB, class TC> inline static void Calc (const TA &a,const TB &b, TC &c) { c = (TC) (a ^ b) ; }	} ;
		class FC_NOT		{ public :	template <class T> inline static void Calc (T a, T &b) { b = !a ; }	} ;
		class FC_not		{ public :	template <class T> inline static void Calc (T a, T &b) { b = ~a ; }	} ;

		class FC_greater		{ public :	template <class TA, class TB, class TC> inline static void Calc (TA a, TB b, TC &c) { c = (TC) (a > b) ; }	} ;
		class FC_greater_equal	{ public :	template <class TA, class TB, class TC> inline static void Calc (TA a, TB b, TC &c) { c = (TC) (a >= b) ; }	} ;
		class FC_equal			{ public :	template <class TA, class TB, class TC> inline static void Calc (TA a, TB b, TC &c) { c = (TC) (a == b) ; }	} ;
		class FC_less		{ public :	template <class TA, class TB, class TC> inline static void Calc (TA a, TB b, TC &c) { c = (TC) (a < b) ; }	} ;
		class FC_less_equal	{ public :	template <class TA, class TB, class TC> inline static void Calc (TA a, TB b, TC &c) { c = (TC) (a <= b) ; }	} ;

		
	} ;

///////////////////////////////////////////////////////////////////////////////
//	CRefBase
///////////////////////////////////////////////////////////////////////////////

	template <class T>
	class CRefBase
	{
		public:

			CRefBase (DWORD dwRefCount = 0) : m_dwRefCount (dwRefCount)
			{
//				Rprintf ("Created Refbase 0x%08p - %d \r\n", this, m_dwRefCount) ;
				ASSERT (IsRefBase (this)) ;		//	ensures, that T is derived from CrefBase :)
			}

			inline static void CutRef ( T * &pref )
			{
//Rprintf ("cutting ref of 0x%08p \r\n", pref) ;
				if (pref)
					pref->DecRef () ;
				pref = NULL ;
			}

			T * Reference ( T * &pref )
			{
//Rprintf ("Referenceing 0x%08p with 0x%08p \r\n", this, pref) ;
				if (this == pref)
					return pref ;

				if (pref)
					pref->DecRef () ;

//				if (m_dwLock)	//	this object is locked (for editing), so it can't be referenced --> make a copy of the object & lock it
//					return (new T (*BACK_CAST (this, T, CRefBase <T>)))->Reference (pref = NULL) ;	//	the constructor must care about deriving from locked objects!

				pref = BACK_CAST (this, T, CRefBase <T>) ;
				IncRef () ;

				return pref ;
			}

			inline BOOL SingleOwner () const { return m_dwRefCount <= 1 ;}
/*
			static T* Empty (T *&pref)
			{
				static T empty ;
				return empty.Reference (pref) ;
//		return CRefBase <T>::Empty		(pref) ;
			}
*/
		protected:
			inline void DecRef ()		
			{
				if (!--m_dwRefCount) 
				{
//					Rprintf ("Deleting Refbase 0x%08p - %d \r\n", this, m_dwRefCount) ;
					delete BACK_CAST (this, T, CRefBase <T>) ;
				}
			}

			inline void IncRef ()		{ m_dwRefCount ++ ; }

			static BOOL IsRefBase (CRefBase <T> *p) { return TRUE ; }
			static BOOL IsRefBase (void *p)			{ return FALSE ; }

			DWORD	m_dwRefCount ;
	} ;

	template <class T> inline void	CutRef		(T *&pref) {		CRefBase <T>::CutRef	(pref) ; }
//	template <class T> inline T *	Edit		(T *&pref) { return CRefBase <T>::Edit		(pref) ; }

#define DEC_EMPTY(TYPE) static TYPE *Empty (TYPE *&pref) { static TYPE empty; return empty.Reference (pref) ; }

/*
	template <class T> inline T *	Empty		(T *&pref)
	{
		return (new T (0))->Reference (pref) ;
//Rprintf ("init emptying\r\n", sizeof (T)) ;
		static T empty ;
//Rprintf ("emptying class of size %d\r\n", sizeof (T)) ;
		T *pret = empty.Reference (pref) ;
//Rprintf ("emptyed class of size %d\r\n", sizeof (T)) ;
		return pret ;
	}
*/
///////////////////////////////////////////////////////////////////////////////
//	CRefBase
///////////////////////////////////////////////////////////////////////////////

	template <class T>
	class CRefLockBase : public CRefBase <T>
	{
		friend class CRefBase <T> ;
		typedef CRefBase <T> t_base ;
		public:

			CRefLockBase (DWORD dwRefCount = 0) : CRefBase<T> (dwRefCount), m_dwLock (0)
			{
				ASSERT (IsRefBase (this)) ;		//	ensures, that T is derived from CrefBase :)
			}

			inline static void CutRef ( T * &pref )
			{
				if (pref)
					pref->DecRef () ;
				pref = NULL ;
			}

			T * Reference ( T * &pref )
			{
				if (this == pref)
					return pref ;

				if (pref)
					pref->DecRef () ;

				if (m_dwLock)	//	this object is locked (for editing), so it can't be referenced --> make a copy of the object & lock it
					return (new T (*BACK_CAST (this, T, CRefBase <T>)))->Reference (pref = NULL) ;	//	the constructor must care about deriving from locked objects!

				pref = BACK_CAST (this, T, CRefBase <T>) ;
				t_base::IncRef () ;

				return pref ;
			}

			T *SoftLock (T *&pref)
			{
				if (this == pref)
					return pref ;

				if (pref)
					Unlock (pref) ;

//				if (!m_dwLock && t_base::m_dwRefCount > 1)
//					return pref = NULL ;

				m_dwLock++ ;
				t_base::IncRef () ;
				pref = BACK_CAST (this, T, CRefLockBase <T>) ;
				return pref ;
			}

			static T *Lock (T *&pref)
			{
				if (pref->m_dwLock || pref->m_dwRefCount <= 1)	//	this class is "lockable" -> lock it
				{
					pref->m_dwLock++ ;
					pref->IncRef () ;
					return pref ;
				}

				return (new T (*pref))->Reference (pref)->Lock (pref) ;
			}

			static BOOL Unlock (T *&pref)
			{
				if (!pref)
					return FALSE ;
				ASSERT (pref->m_dwLock >= 1) ;
				BOOL bRet = !pref->m_dwLock || --pref->m_dwLock ;
				pref->DecRef () ;
				pref = NULL ;
				return bRet ;
			}

			inline BOOL IsLocked () const { return m_dwLock > 0 ; }
		protected:

			DWORD	m_dwLock ;
	} ;

	template <class T> inline T *	Lock		(T *&pref) { return CRefLockBase <T>::Lock		(pref) ; }
	template <class T> inline BOOL	Unlock		(T *&pref) { return CRefLockBase <T>::Unlock	(pref) ; }

///////////////////////////////////////////////////////////////////////////////
//	CDataRef
///////////////////////////////////////////////////////////////////////////////

	enum SetDataFlag
	{
		sdf_modeReference = 0,
		sdf_modeAttach = 1,
		sdf_modeCopy = 2
	} ;

class TClust ;
	template <class T>
	class CDataRef : public CRefLockBase < CDataRef <T> >
	{
		friend class TClust ; 
		friend class TClust1 ;
		friend class CRefBase <T> ;
		friend class CRefLockBase <T> ;
		typedef CRefLockBase < CDataRef <T> > t_base ;
	public:
		DEC_EMPTY (CDataRef <T> )

		CDataRef (T *pData = NULL, DWORD dwSize = 0, SetDataFlag sdf = sdf_modeReference, void *pOwner = NULL) : t_base (1), m_pData (0), m_dwSize (0), m_pOwner (NULL)
		{
			SO_SetData (dwSize, pData, sdf, pOwner) ;
		}

		~CDataRef ()
		{
			ASSERT (!m_dwLock) ;	//	if this dataref MUST NOT be locked on destruction - this must be ensured by the programmer
	
			if (!m_pOwner)
			{
//				::Rprintf ("freeing 0x%08p ...", m_pData) ;
				free (m_pData) ;
//				::Rprintf (" done! \r\n") ;

			}
		}

		static CDataRef<T> *Create (DWORD dwSize, CDataRef<T> *&pref, T *pData = NULL, SetDataFlag sdf = sdf_modeReference, void *pOwner = NULL)
		{
			if (!pData)
				return (new CDataRef<T> (0, (T *) malloc (dwSize * sizeof (T)), dwSize, sdf_modeAttach))->Reference (pref) ;
			return (new CDataRef<T> (0, pData, dwSize, sdf, pOwner))->Reference (pref) ;
		}

		static CDataRef<T> *ReCreate (DWORD dwSize, CDataRef<T> *&pref, T *pData = NULL, SetDataFlag sdf = sdf_modeReference, void *pOwner = NULL)
		{
			if (!pref || pref->m_dwRefCount > 1)
				return Create (dwSize, pref, pData, sdf, pOwner) ;

			if (pData)
				pref->SO_SetData (dwSize, pData, sdf, pOwner) ;
			else
				pref->SO_SetSize (dwSize) ;
			return pref ;				
		}

		static T *Detach (CDataRef <T> *&pref)
		{
			if (!pref->SingleOwner ())
				( new CDataRef <T> (*pref) )->Reference  (pref) ;

			pref->m_dwSize = 0 ;

			T *pRet = pref->m_pData ;
			pref->m_pData = NULL ;

			Empty (pref) ;		//	now the old object is destroyed and the refernce pointer is set to an (the) emtpy object
			return pRet ;
		}

		inline T &GetItem (DWORD dwIdx) { ASSERT (dwIdx < m_dwSize) ; 	return m_pData [dwIdx] ; }

		inline DWORD GetSize ()	const	{ return m_dwSize ; }
		inline DWORD GetMemSize ()	const	{ return m_dwSize * sizeof (T) ; }
		inline T *GetData () { return m_pData ; }

		template <class U>
		CDataRef <U> * &Convert (CDataRef <U> *&pref)
		{
			pref->Create (m_dwSize, pref) ;

			U *pCurU = pref->GetData () ;
			T *pCurT = GetData () ;

			DWORD i ;
			for (i = m_dwSize; i; i--)
				*pCurU++ = (U) *pCurT++ ;

			return pref ;
		}

		CDataRef (DWORD dwRefCount, T *pData, DWORD dwSize, SetDataFlag sdf = sdf_modeReference, void *pOwner = NULL) : t_base (dwRefCount), m_pData (NULL), m_dwSize (0), m_pOwner (NULL)
		{
			SO_SetData (dwSize, pData, sdf, pOwner) ;
		}

		CDataRef (const CDataRef <T> &ref)
		{
			SO_SetData (ref.m_dwSize, ref.m_pData, sdf_modeCopy) ;
		}

		T *CreateCopyDetach ()
		{		//	copies the referenced data (such, that it can be used by all objects referencing this class) and returns the pointer to the original data
			T *pRet = m_pData ;

			m_pData = (T *) malloc (GetMemSize ()) ;
			memcpy (m_pData, pRet, GetMemSize ()) ;
			m_pOwner = NULL ;

			return pRet ;
		}

		DWORD GetRefCount () { 	return CRefBase < CDataRef <T> >::m_dwRefCount ; }

		void *GetOwner () const { return m_pOwner ; }

	protected:

		BOOL SO_SetSize (DWORD dwSize)
		{
			ASSERT (SingleOwner ()) ;

			if (dwSize > m_dwSize || dwSize < (m_dwSize >> 1))
				SO_ForceSize (dwSize) ;
			return TRUE ;
		}

		void SO_ForceSize (DWORD dwSize)
		{
			ASSERT (SingleOwner ()) ;

			m_dwSize = dwSize ;
			if (m_pOwner)
			{
				m_pData = (T *) malloc (GetMemSize ()) ;
				m_pOwner = NULL ;
			}
			else
				m_pData = (T *) realloc (m_pData, GetMemSize ()) ;
		}

		void SO_SetData (DWORD dwSize, T *pData = NULL, SetDataFlag sdf = sdf_modeReference, void *pOwner = NULL)
		{			//	setzt single Owner vorraus
			ASSERT (SingleOwner ()) ;

			if (!m_pOwner)
				free (m_pData) ;


			m_pOwner = NULL ;
			m_dwSize = dwSize ;
			if (sdf == sdf_modeCopy)
			{
				m_pData = (T *) malloc (GetMemSize ()) ;
				memcpy (m_pData, pData, GetMemSize ()) ;
			}
			else
			{
				m_pData = pData ;
				if (sdf == sdf_modeReference)
					m_pOwner = pOwner ;
			}
		}

		T		*m_pData ;
		DWORD	m_dwSize ;
		void 	*m_pOwner ;
	} ;

	template <class T> inline T *Detach (CDataRef <T> *&pref) { return CDataRef<T>::Detach (pref) ; }

/*
///////////////////////////////////////////////////////////////////////////////
//	CEditBase 
///////////////////////////////////////////////////////////////////////////////


	template <class T>
	class CEditBase
	{
		public:
			CEditBase (const CEditBase &editbase) : m_pEditItem (editbase.m_pEditItem)
			{
				ASSERT (m_pEditItem->IsLocked ()) ;
				Lock (m_pEditItem) ;
			}

			CEditBase (T * &pEditItem)
			{
				m_pEditItem = Lock (pEditItem) ;
			}

			~CEditBase ()
			{
				if (m_pEditItem)
					m_pEditItem->Unlock () ;
			}

			CEditBase<T> & operator = (CEditBase<T> &editbase)
			{
				if (m_pEditItem)
					m_pEditItem->Unlock () ;
				m_pEditItem = editbase.m_pEditItem ;
				Lock (m_pEditItem) ;
			}

			BOOL StopEdit ()
			{
				BOOL bRet = m_pEditItem->Unlock () ;
				m_pEditItem = NULL ;
				return bRet ;
			}

		protected:

			T *GetEditItem () { return m_pEditItem ; }

			T *m_pEditItem ;
	} ;
*/

///////////////////////////////////////////////////////////////////////////////
//	CRefRef
///////////////////////////////////////////////////////////////////////////////

	template <class T>
	class CRefRef
	{
	public:
		CRefRef (T *pItem) { pItem->Reference (m_pItem = 0) ; }
		CRefRef (CRefRef <T> *pRefRef) { pRefRef->m_pItem->Reference (m_pItem = 0) ; }

		~CRefRef () { CutRef (m_pItem) ; }
		T *GetItem () { return m_pItem ; }

	protected:
		T *m_pItem ;
	} ;

///////////////////////////////////////////////////////////////////////////////
//	CSimpleArray
///////////////////////////////////////////////////////////////////////////////

	template <class T>
	class CSimpleArray
	{
		public:

		CSimpleArray (DWORD dwSize) : m_pData (new T [dwSize]) {}
		CSimpleArray () : m_pData (NULL) {}

		~CSimpleArray () { delete [] m_pData ; }

		void Create (DWORD dwSize)
		{
			delete [] m_pData ;
			m_pData = new T [dwSize] ;
#ifdef _DEBUG 
			m_dwSize = dwSize ;
#endif
		}

		inline T &operator [] (DWORD d)
		{
			ASSERT (d < m_dwSize) ;
			return m_pData[d] ;
		}

		inline const T &operator [] (DWORD d) const
		{
			ASSERT (d < m_dwSize) ;
			return m_pData[d] ;
		}


		protected:
		T *m_pData ;
#ifdef _DEBUG
		DWORD m_dwSize ;
#endif

	} ;



