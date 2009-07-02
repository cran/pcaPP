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

#include "IRef.h"
#include "ElOp.h"
#include "stdio.h"
#include "stdlib.h"
#include "IMat.h"
#include "IDim.h"

//	class ISub ;
//	class ISubR ;

//	template <class T> class IVec ;
//	template <class T> class IMat ;

/*
	template <class T>
	BOOL rbind (IMat <T> &res, const char * szInput, ... )
	{
		
	}

	template <class T>
	BOOL cbind (IMat <T> &res, const char * szInput, ... )
	{
		
	}
*/

	template <class PAR1, class PAR2>
	void sortmatrow (PAR1 &matin, PAR2 &matout, IVec<DWORD> &idx)
	{
//		matout.Create (idx.size (), matin.ncol ()) ;
//		typename PAR2::m_edittype retedit (matout .Edit ()) ;

		DWORD r, c ;
		for (r = idx.size () - 1; r != (DWORD) -1; r--)
		{
			DWORD dwCurR = idx(r) ;
			for (c = matin.ncol () - 1; c != (DWORD) -1; c--)
				ret(r, c) = matin (dwCurR, c) ;
		}
//		retedit.StopEdit () ;
	}

	template <class PAR1, class PAR2>
	void sortmatcol (PAR1 &matin, PAR2 &matout, IVec<DWORD> &idx)
	{
//		matout.Create (matin.nrow (), idx.size ()) ;
//		typename PAR2::m_edittype retedit (matout .Edit ()) ;

		DWORD r, c ;
		for (c = idx.size () - 1; c != (DWORD) -1; c--)
		{
			DWORD dwCurC = idx(c) ;
			for (r = matin.nrow () - 1; r != (DWORD) -1; r--)
				ret(r, c) = matin (r, dwCurC) ;
		}
//		retedit.StopEdit () ;
	}

	template <class PAR1, class PAR2>
	void sortvec (PAR1 &vecin, PAR2 &vecout, IVec <DWORD> &idx)
	{
//		vecout.Create (idx.size ()) ;
//		typename PAR2::m_edittype retedit (vecout.Edit ()) ;

		DWORD v ;
		for (v = idx.size () - 1; v != (DWORD) -1; v--)
			ret (v) = vecin (idx(v)) ;

//		retedit.StopEdit () ;
	}

	template <class RET>
	class DimOrder
	{

		template <class PAR>
		IMat<RET> sortmatrow (PAR &mat, IVec<DWORD> &idx)
		{
			IMat<RET> ret (idx.size (), mat.ncol ()) ;
			::sortmatrow (mat, ret, idx) ;
			return ret ;
		}

		template <class PAR>
		IMat<RET> sortmatcol (PAR &mat, IVec<DWORD> &idx)
		{
			IMat<RET> ret (mat.nrow (), idx.size ()) ;
			::sortmatcol (mat, ret, idx) ;
			return ret ;
		}

		template <class PAR>
		IVec<RET> sortvec (PAR &vec, IVec <DWORD> &idx)
		{
			IVec<RET> ret (idx.size()) ;
			::sortvec (vec, ret, idx) ;
			return ret ;
		}
	} ;

