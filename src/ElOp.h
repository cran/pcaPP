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


	template <class T> class IMat ;
	template <class T> class IVec ;
//	template <class T> class IMatEdit ;
//	template <class T> class IVecEdit ;
	template <class T, class D> class ITens_ ;
//	template <class T, class D> class ITensEdit_ ;
	template <class T, class D> class ITensConst ;

	inline DWORD CheckIdx (DWORD dwIdx, DWORD dwSize) { return (dwIdx < dwSize) ? dwIdx : dwIdx % dwSize ; }

	template <class F, class RET>
	class FC_ElOp
	{
		protected: 
		public :
		template <class PAR1, class PAR2>
		static IMat<RET> OpMM (const PAR1 &mat1, const PAR2 &mat2, DWORD dwRetRow = 0, DWORD dwRetCol = 0)	//	elementary operation (2matrices)
		{
			if (!dwRetRow)
				dwRetRow = max (mat1.nrow (), mat2.nrow ()) ;
			if (!dwRetCol)
				dwRetCol = max (mat1.ncol (), mat2.ncol ()) ;

			DWORD	r1 = CheckIdx (dwRetRow - 1, mat1.nrow ()),
					r2 = CheckIdx (dwRetRow - 1, mat2.nrow ()),
					c1 = CheckIdx (dwRetCol - 1, mat1.ncol ()),
					c2 = CheckIdx (dwRetCol - 1, mat2.ncol ()) ;

			const DWORD	nrow1m1 = mat1.nrow () - 1, nrow2m1 = mat2.nrow () - 1, 
						ncol1m1 = mat1.ncol () - 1, ncol2m1 = mat2.ncol () - 1 ;

			IMat<RET> ret (dwRetRow, dwRetCol) ;
//			IMatEdit<RET> retedit (ret.Edit ()) ;

			DWORD r, c ;
			for (r = dwRetRow - 1; r != (DWORD) -1; r--)
			{
				for (c = dwRetCol - 1; c != (DWORD) -1; c--)
				{
					F::Calc (mat1 (r1, c1), mat2 (r2, c2), ret (r, c)) ;
					if (--c1 == (DWORD) -1)	c1 = ncol1m1 ;
					if (--c2 == (DWORD) -1)	c2 = ncol2m1 ;
				}
				if (--r1 == (DWORD) -1)	r1 = nrow1m1 ;
				if (--r2 == (DWORD) -1)	r2 = nrow2m1 ;
			}

//			retedit.StopEdit () ;
			return ret ;
		}

		template <class PAR>
		static IMat<RET> &OpM (PAR &mat)
		{
			IMat<RET> ret (mat.nrow (), mat.ncol ()) ;
//			IMatEdit<RET> retedit = mat.Edit () ;
			
			const DWORD	ncolm1 = mat.ncol () - 1 ;

			DWORD r, c ;
			for (r = mat.nrow () - 1; r != (DWORD) -1; r--)
				for (c = ncolm1; c != (DWORD) -1; c--)
					F::Calc (ret (r, c), mat (r, c)) ;
			return ret ;
		}

		template <class PAR>
		static IVec<RET> &OpV (PAR &vec)
		{
			IVec<RET> ret (vec.nrow (), vec.ncol ()) ;
//			IVecEdit<RET> retedit = vec.Edit () ;

			DWORD v = vec.size () - 1;
			for (; v != (DWORD) -1; v--)
				F::Calc (ret (v), vec (v)) ;
			return vec ;
		}

		template <class PAR1, class PAR2>
		static IMat<RET> OpMV (const PAR1 &mat, const PAR2 &vec)
		{
			IMat<RET> ret (mat.nrow (), mat.ncol ()) ;
//			IMatEdit<RET> retedit = ret.Edit () ;

			DWORD dwVecIdx = (mat.nrow () * mat.ncol () - 1) % vec.size () ;

			const DWORD	dwMatnRowm1 = mat.nrow () - 1, 
						dwVecSizem1 = vec.size () - 1 ;

			DWORD r, c ;
			for (c = mat.ncol () - 1; c != (DWORD) -1; c--)
				for (r = dwMatnRowm1; r != (DWORD) -1; r--)
				{
					F::Calc (mat (r, c), vec (dwVecIdx--), ret (r, c)) ;
					if (dwVecIdx == (DWORD) -1)
						dwVecIdx = dwVecSizem1 ;
				}

//			retedit.StopEdit () ;
			return ret ;
		}

		template <class PAR1, class PAR2>
		static const IMat<RET> OpMV_col (const PAR1 &mat, const PAR2 &vec)
		{
			IMat<RET> ret (mat.nrow (), mat.ncol ()) ;
			return OpMV_col (mat, vec, ret) ;
		}


		template <class T, class U, class V>
		static const IMat<T> &OpMV_col (const IMat<T> &mat, const IVec <U> &vec, const IMat<V> &ret)
		{
			DWORD dwVecIdx = (mat.nrow () * mat.ncol () - 1) % vec.size () ;

			const DWORD	dwMatnRowm1 = mat.nrow () - 1, 
						dwVecSizem1 = vec.size () - 1 ;

			DWORD r, c ;
			for (c = mat.ncol () - 1; c != (DWORD) -1; c--)
				for (r = dwMatnRowm1; r != (DWORD) -1; r--)
				{
					double &d1 = mat (r, c) ;
					double &d2 = vec (dwVecIdx) ;
					double &dret = ret (r, c) ;

					F::Calc (d1, d2, dret) ;

					dwVecIdx-- ;

					if (dwVecIdx == (DWORD) -1)
						dwVecIdx = dwVecSizem1 ;
				}


			return ret ;
		}

		template <class PAR1, class PAR2>
		static const IMat<RET> OpMV_row (const PAR1 &mat, const PAR2 &vec)
		{
			IMat<RET> ret (mat.nrow (), mat.ncol ()) ;
			return OpMV_row (mat, vec, ret) ;
		}


		template <class T, class U, class V>
		static const IMat<T> &OpMV_row (const IMat<T> &mat, const IVec <U> &vec, const IMat<V> &ret)
		{
			DWORD dwVecIdx = (mat.nrow () * mat.ncol () - 1) % vec.size () ;

			const DWORD	dwMatnColm1 = mat.ncol () - 1, 
						dwVecSizem1 = vec.size () - 1 ;

			DWORD r, c ;
			for (r = mat.nrow () - 1; r != (DWORD) -1; r--)
				for (c = dwMatnColm1; c != (DWORD) -1; c--)
				{
					double &d1 = mat (r, c) ;
					double &d2 = vec (dwVecIdx) ;
					double &dret = ret (r, c) ;

					F::Calc (d1, d2, dret) ;

					dwVecIdx-- ;


					if (dwVecIdx == (DWORD) -1)
						dwVecIdx = dwVecSizem1 ;
				}


			return ret ;
		}

		template <class PAR1, class PAR2>
		static IMat<RET> OpVM (const PAR1 &vec, const PAR2 &mat)
		{
			IMat<RET> ret (mat.nrow (), mat.ncol ()) ;
//			IMatEdit<RET> retedit (ret.edit ()) ;

			DWORD v = (mat.nrow () * mat.ncol ()) % vec.size () ;

			const DWORD	dwMatnRowm1 = mat.nrow () - 1, 
						dwVecSizem1 = vec.size () - 1 ;

			DWORD r, c ;
			for (c = mat.ncol () - 1; c != (DWORD) -1; c--)
				for (r = dwMatnRowm1; r != (DWORD) -1; r--)
				{
					F::Calc (vec (v--), mat (r, c), ret (r, c)) ;
					if (v == (DWORD) -1)	v = dwVecSizem1 ;
				}

//			retedit.StopEdit () ;
			return ret ;
		}

		template <class PAR1, class PAR2>
		static IMat<RET> OpVM_row (const PAR1 &vec, const PAR2 &mat)
		{
			IMat<RET> ret (mat.nrow (), mat.ncol ()) ;
//			IMatEdit<RET> retedit ( ret.edit ()) ;

			DWORD v = (mat.nrow () * mat.ncol ()) % vec.size () ;

			const DWORD	dwMatnColm1 = mat.ncol () - 1, 
					dwVecSizem1 = vec.size () - 1 ;

			DWORD r, c ;
			for (r = mat.nrow () - 1; r != (DWORD) -1; r--)
				for (c = dwMatnColm1; c != (DWORD) -1; c--)
				{
					F::Calc (vec (v--), mat (r, c), ret (r, c)) ;
					if (v == (DWORD) -1)	v = dwVecSizem1 ;
				}

///			retedit.StopEdit () ;
			return ret ;
		}

		template <class PAR1, class PAR2>
		static IVec<RET> OpVV (const PAR1 &vec1, const PAR2 &vec2, DWORD dwRetSize = 0)
		{
			if (!dwRetSize)
				dwRetSize = max (vec1.size (), vec2.size ()) ;

			IVec<RET> ret (dwRetSize) ;
			return OpVV (vec1, vec2,ret) ;
//			IVecEdit<RET> retedit = ret.edit () ;

/*			DWORD	dwIdx1 = CheckIdx (dwRetSize - 1, vec1.size ()),
					dwIdx2 = CheckIdx (dwRetSize - 1, vec2.size ()) ;
		
			DWORD i ;
			for (i = dwRetSize - 1; i != (DWORD) -1; i--)
			{
				F::Calc (vec1 (dwIdx1--), vec2 (dwIdx2--), ret (i)) ;
				if (dwIdx1 == (DWORD) -1)
					dwIdx1 = vec1.size () - 1 ;
				if (dwIdx2 == (DWORD) -1)
					dwIdx2 = vec2.size () - 1 ;
			}

//			retedit.StopEdit () ;
			return ret ;
*/		}

		template <class PAR1, class PAR2>
		static IVec<RET> &OpVV (const IVec <PAR1> &vec1, const IVec <PAR2> &vec2, IVec<RET> &ret)
		{

			DWORD dwRetSize = max (vec1.size (), vec2.size ()) ;
			ret.Reshape (dwRetSize) ;

			DWORD	dwIdx1 = CheckIdx (dwRetSize - 1, vec1.size ()),
					dwIdx2 = CheckIdx (dwRetSize - 1, vec2.size ()) ;
		
			DWORD i ;
			for (i = dwRetSize - 1; i != (DWORD) -1; i--)
			{
				F::Calc (vec1 (dwIdx1--), vec2 (dwIdx2--), ret (i)) ;
				if (dwIdx1 == (DWORD) -1)
					dwIdx1 = vec1.size () - 1 ;
				if (dwIdx2 == (DWORD) -1)
					dwIdx2 = vec2.size () - 1 ;
			}

			return ret ;
		}

		template <class PAR1, class U>
		static IMat<RET> OpME (const PAR1 &mat, U val)
		{
			IMat<RET> ret (mat.ncol (), mat.nrow()) ;
//			IMatEdit<RET> retedit = ret.Edit () ;

			DWORD r, c ;
			for (r = ret.nrow () - 1; r != (DWORD) -1; r--)
				for (c = ret.ncol () - 1; c != (DWORD) -1; c--)
					 F::Calc (mat (r, c), val, ret (r, c)) ;
//			retedit.StopEdit () ;
			return ret ;
		}

		template <class PAR1, class U>
		static const IVec<RET> OpVE (const PAR1 &vec, U val)
		{
			IVec<RET> ret (vec.size ()) ;
			return OpVE (vec, val, ret) ;
		}

		template <class T, class U, class V>
		static const IVec<V> &OpVE (const IVec <T> &vec, const U &val, const IVec<V> &ret)
		{
			DWORD i ;
			for (i = vec.size () - 1; i != (DWORD) -1; i--)
				F::Calc (vec (i),
				
				val,
				
				ret (i)) ;

			return ret ;
		}

		template <class T, class U, class V>
		static V *OpPP (const T *pS1, const DWORD dwSize1, const U *pS2, DWORD dwSize2, V *pD, DWORD dwSizeD)
		{
			DWORD	c1 = (dwSizeD -  1) % dwSize1,
					c2 = (dwSizeD -  1) % dwSize2,
					c ;

			for (c = dwSizeD - 1; c != (DWORD) -1; c--)
			{
				F::Calc (pS1[c1--], pS2[c2--], pD[c]) ;
				if(c1 == (DWORD) -1) c1 = dwSize1 - 1 ;
				if(c2 == (DWORD) -1) c2 = dwSize2 - 1 ;
			}
			return pD ;
		}

		template <class T, class U>
		static RET *OpPE (const T *	pS, const DWORD dwSize, const U val, RET *pD)
		{
			DWORD c ;
			for (c = dwSize - 1; c != (DWORD) -1; c--)
				F::Calc (pS[c], val, pD[c]) ;
			return pD ;
		}


		
/*		template <class T, class U, class V, class D>
		static ITens_<V, D> 
		&OpTE 
		(const T &tens, 
		const U &val, 
		ITens_<V, D> &ret)
		{
			OpTE (tens, val, ret).StopEdit () ;
			return ret ;
		}
*/
		template <class T, class U, class V, class D>
		static const ITens_<V, D> &OpTE (const T &tens, const U &val, const ITens_<V, D> &ret)
		{
			typename ITens_<V, D>::t_Iter	dim3 = ret.dim (),
								iter3 = dim3 - (DWORD) 1 ;
			typename T::t_Iter	dim1 = tens.dim (),
								iter1 = iter3 % dim1 ;
			do
			{
				F::Calc (tens(iter1), val, ret(iter3)) ;
			} while (iter3.mm (dim3)) ;

			return ret ;
		}



	} ;

	template <class F>
	class FC_ElOpAs
	{
		public:
		template <class PAR1, class PAR2>
		static PAR1 &OpMM (PAR1 &mat1, const PAR2 &mat2)	//	elementary operation
		{
//			typename PAR1::m_edittype mat1edit = mat1.Edit () ;

			DWORD	r2 = CheckIdx (mat1.nrow () - 1, mat2.nrow ()),
					c2 = CheckIdx (mat1.ncol () - 1, mat2.ncol ()) ;

			const DWORD	ncol1m1 = mat1.ncol () - 1,
						ncol2m1 = mat2.ncol () - 1,
						nrow2m1 = mat2.nrow () - 1 ;

			DWORD r, c ;
			for (r = mat1.nrow () - 1; r != (DWORD) -1; r--)
			{
				for (c = ncol1m1; c != (DWORD) -1; c--)
				{
					F::Calc (mat1 (r, c), mat2 (r2, c2), mat1 (r, c)) ;
					if (--c2 == (DWORD) -1)	c2 = ncol2m1 ;
				}
				if (--r2 == (DWORD) -1)	r2 = nrow2m1 ;
			}

//			mat1edit.StopEdit () ;
			return mat1 ;
		}

		template <class PAR>
		static PAR &OpM (PAR &mat)
		{
//			typename PAR::m_edittype matedit = mat.Edit () ;
			
			const DWORD	ncolm1 = mat.ncol () - 1 ;

			DWORD r, c ;
			for (r = mat.nrow () - 1; r != (DWORD) -1; r--)
				for (c = ncolm1; c != (DWORD) -1; c--)
					F::Calc (mat (r, c), mat (r, c)) ;
			return mat ;
		}

		template <class PAR>
		static PAR &OpV (PAR &vec)
		{
//			typename PAR::m_edittype vecedit = vec.Edit () ;
			
			DWORD v = vec.size () - 1;
			for (; v != (DWORD) -1; v--)
				F::Calc (vec (v), vec (v)) ;
			return vec ;
		}

		template <class PAR1, class PAR2>
		static const IMat <PAR1> &OpMV (const IMat <PAR1> &mat, const IVec <PAR2> &vec)
		{
//			typename PAR1::m_edittype matedit = mat.edit () ;

			DWORD dwVecIdx = (mat.nrow () * mat.ncol ()) % vec.size () ;

			const DWORD	dwMatnRowm1 = mat.nrow () - 1, 
						dwVecSizem1 = vec.size () - 1 ;

			DWORD c, r ;
			for (c = mat.ncol () - 1; c != (DWORD) -1; c--)
				for (r = dwMatnRowm1; r != (DWORD) -1; r--)
				{
					F::Calc (mat (r, c), vec (dwVecIdx--), mat (r, c)) ;
					if (dwVecIdx == (DWORD) -1)	dwVecIdx = dwVecSizem1 ;
				}

//			matedit.StopEdit () ;
			return mat ;
		}


		template <class PAR1, class PAR2>
		static PAR1 &OpMV_row (PAR1 &mat, const PAR2 &vec)
		{
//			typename PAR1::m_edittype matedit = mat.edit () ;

			DWORD dwVecIdx = (mat.nrow () * mat.ncol () - 1) % vec.size () ;

			const DWORD	dwMatnColm1 = mat.ncol () - 1, 
						dwVecSizem1 = vec.size () - 1 ;

			DWORD r, c ;
			for (r = mat.nrow () - 1; r != (DWORD) -1; r--)
				for (c = dwMatnColm1; c != (DWORD) -1; c--)
				{
					F::Calc (mat (r, c), vec (dwVecIdx), mat (r, c)) ;
					if (--dwVecIdx == (DWORD) -1)	dwVecIdx = dwVecSizem1 ;
				}

//			matedit.StopEdit () ;
			return mat ;
		}

		template <class PAR1, class U>
		static PAR1 &OpME (PAR1 &mat, const U val)	//	elementary operation
		{
//			typename PAR1::m_edittype matedit = mat.Edit () ;

			int nCol1m1 = mat.ncol () - 1 ;

			DWORD r, c ;
			for (r = mat.nrow () - 1; r != (DWORD) -1; r--)
				for (c = nCol1m1; c != (DWORD) -1; c--)
					F::Calc (mat (r, c), val, mat (r, c)) ;

//			matedit.StopEdit () ;
			return mat ;
		}

		template <class T, class U>
		static const IVec<T> &OpVV (const IVec <T> &vec1, const IVec <U> &vec2)
		{
//			typename PAR1::m_edittype vecedit =7 vec1.edit () ;

			DWORD	i = vec1.size () -  1 ;
			DWORD	dwIdx2 = CheckIdx (i, vec2.size ()) ;
		
			for (i = i; i != (DWORD) -1; i--)
			{
				F::Calc (vec1 (i), vec2 (dwIdx2--), vec1 (i)) ;
				if (dwIdx2 == (DWORD) -1)
					dwIdx2 = vec2.size () - 1 ;
			}

//			vecedit.StopEdit () ;
			return vec1 ;
		}

		template <class PAR1, class U>
		static PAR1 &OpVE (PAR1 &vec, const U val)	//	elementary operation
		{
//			typename PAR1::m_edittype vecedit = vec.Edit () ;

			DWORD i ;
			for (i = vec.size () - 1; i != (DWORD) -1; i--)
				F::Calc (vec (i), val, vec (i)) ;

//			vecedit.StopEdit () ;
			return vec ;
		}

		template <class T, class U>
		static T *OpPP (const T *pS1, const DWORD dwSize1, const U *pS2, DWORD dwSize2)
		{
			DWORD	c = dwSize1 - 1,
					c2 = c % dwSize2 ;

			for (; c != (DWORD) -1; c--)
			{
				F::Calc (pS1[c], pS2[c2--], pS1[c]) ;
				if(c2 == (DWORD) -1) c2 = dwSize2 - 1 ;
			}

			return pS1 ;
		}

		template <class T, class U>
		static T *OpPE (const T *pS, const DWORD dwSize, const U val)
		{
			DWORD c ;
			for (c = dwSize - 1; c != (DWORD) -1; c--)
				F::Calc (pS[c], val, pS[c]) ;
			return pS ;
		}

		template <class T, class U, class D>
		static ITens_<T, D> &OpTE (ITens_<T, D> &tens, const U &val)
		{
			OpTE (tens, val) ;
			return tens ;
		}

		template <class T, class U, class D>
		static const ITens_<T, D> &OpTE (const ITens_<T, D> &tens, const U &val)
		{
			typename ITens_<T, D>::t_Iter	dim = tens.dim (),
											iter = dim - (DWORD) 1 ;
			do
			{
				typename ITens_<T, D>::m_datatype &cur = tens(iter) ;
				F::Calc (cur, val, cur) ;
			} while (iter.mm (dim)) ;

			return tens ;
		}
	} ;

//	template <class F>
	class FC_As
	{
		public:
		template <class T, class PAR2>
		static const IMat <T> &OpMM (const IMat <T> &mat1, const PAR2 &mat2)	//	elementary operation
		{
//			typename PAR1::m_edittype mat1edit = mat1.Edit () ;

			const DWORD	ncol1m1 = mat1.ncol () - 1,
						ncol2m1 = mat2.ncol () - 1,
						nrow2m1 = mat2.nrow () - 1 ;

			DWORD	r2 = CheckIdx (mat1.nrow () - 1, mat2.nrow ()),
					c2 = CheckIdx (ncol1m1, mat2.ncol ()) ;

			DWORD r, c ;
			for (r = mat1.nrow () - 1; r != (DWORD) -1; r--)
			{
				for (c = ncol1m1; c != (DWORD) -1; c--)
				{
					mat1 (r, c) =  mat2 (r2, c2) ;
					if (--c2 == (DWORD) -1)	c2 = ncol2m1 ;
				}
				if (--r2 == (DWORD) -1)	r2 = nrow2m1 ;
			}

//			mat1edit.StopEdit () ;
			return mat1 ;
		}

/*		template <class T, class D, class U, class E>
		static ITens_ <T, D> &OpTT (ITens_ <T, D> &t1, const ITensConst <U, E> &t2)
		{
			OpTT (t1.Edit (), t2) ;
			return t1 ;
		}
*/

		template <class T, class D, class U, class E>
		static const ITens_ <T, D> &OpTT (const ITens_ <T, D> &t1, const ITens_ <U, E> &t2)
		{
			typename ITens_ <T, D>::t_Iter
				itdim1 = t1.dim (),
				it1 = itdim1 - (DWORD) 1 ;

			typename ITens_ <U, E>::t_Iter
				itdim2 = t2.dim (),
				it2 = it1 % itdim2 ;

			do
			{
				t1 (it1) = (T) t2 (it2) ;
				it2.mm (itdim2) ;
			}
			while (it1.mm (itdim1)) ;

			return t1 ;
		}
/*
		template <class T, class D, class U>
		static ITens_ <T, D> &OpTR (ITens_ <T, D> &tens, const U *r)
		{
			OpTR (tens.Edit (), r).StopEdit () ;
			return tens ;
		}
*/

		template <class T, class D, class U>
		static const ITens_ <T, D> &OpTR (const ITens_ <T, D> &tens, const U *r)
		{
			ASSERT (tens.size ()) ;
			typename ITens_ <T, D>::t_Iter	itdim = tens.dim (),
												it = itdim - (DWORD) 1 ;

			r += tens.size () ;

			do
			{
				tens (it) = (T) *--r ;
			}
			while (it.mm (itdim)) ;

			return tens ;
		}

		template <class PAR1, class T>
		static T *OpRT (T *r, const PAR1 &t)
		{
			ASSERT (t.size ()) ;
			typename PAR1::t_Iter	itdim = t.dim (),
									it = itdim - (DWORD) 1 ;

			r += t.size () ;

			do
			{
				*--r = (T) t (it) ;
			}
			while (it.mm (itdim)) ;

			return r ;
		}

		template <class PAR, class T>
		static PAR &OpTE (PAR &t, const T &val)
		{
//			ASSERT (t.size ()) ;	//	this function already handles 
			typename PAR::t_Iter	it, itdim = t.dim () ;
			if (!itdim.Dim2Iter (it))
				return t ;

			while (it.mm (itdim))
			{
				t (it) = (typename PAR::m_datatype) val ;
			}
			

			return t ;
		}

		template <class T, class PAR>
		static const IVec <T> &OpVV (const IVec <T> &vec1, const PAR &vec2)
		{
			ASSERT (vec1.size ()) ;
			DWORD	dwIdx1 = vec1.size () - 1 ,
					dwIdx2 = CheckIdx (dwIdx1, vec2.size ()) ;
					
			for (; dwIdx1 != (DWORD) -1; dwIdx1--)
			{
				vec1(dwIdx1) = vec2 (dwIdx2--) ;
				if (dwIdx2 == (DWORD) -1)
					dwIdx2 = vec2.size () - 1 ;
			}

			return vec1 ;
		}

	} ;


	class FC_mult
	{
		template <class TA, class TB, class TTC>
		BOOL matmultmat (const TA &a, const TB &b, const IMat<TTC> &res)
		{
			ASSERT (!a.IsLocked () && !b.IsLocked ()) ;
			ASSERT (a.ncol () == b.nrow ()) ;

			if (a.ncol () != b.nrow ())
				return FALSE ;		//	make an assertion!!! - check assertions!

	//		res.Create (a.nrow (), b.ncol ()) ;
	//		IMatEdit<TTC> resedit = res.Edit () ;
	//		resedit.Reset () ;

			DWORD dwBncolm1 = b.ncol () - 1 ;
			DWORD dwAncolm1 = a.ncol () - 1 ;

			DWORD i, j, h ;
			for (i = a.nrow () - 1; i != (DWORD) -1; i--)
				for (j = dwBncolm1; j != (DWORD) -1; j--)
				{
					TTC &cur = res (i, j) = 0 ;
					for (h = dwAncolm1; h != (DWORD) -1; h--)
						 cur += (TTC) (a (i, h) * b(h, j)) ;
				}

			return TRUE ;
		}

		template <class TA, class TB, class TTC>
		BOOL vecmultmat (TA &a, TB &b, IMat<TTC> &res)
		{
			ASSERT (!a.IsLocked () && !b.IsLocked ()) ;
			ASSERT (a.size () == b.nrow ()) ;

			if (a.size () != b.nrow ())
				return FALSE ;		//	make an assertion!!! - check assertions!

			res.Create (1, b.ncol ()) ;

//			IMatEdit<TTC> resedit = res.Edit () ;

			res.Reset () ;

			DWORD dwAncolm1 = a.size () - 1 ;

			DWORD j, h ;
			for (j = b.ncol () - 1 ; j != (DWORD) -1; j--)
				for (h = dwAncolm1; h != (DWORD) -1; h--)
					res (0, j) += (TTC) (a (h) * b(h, j)) ;

			return TRUE ;
		}

		template <class TA, class TB, class TTC>
		BOOL matmultvec (TA &a, TB &b, IMat<TTC> &res)
		{
			ASSERT (!a.IsLocked () && !b.IsLocked ()) ;
			ASSERT (a.ncol () == b.size ()) ;

			if (a.ncol () != b.size ())
				return FALSE ;		//	make an assertion!!! - check assertions!

			res.Create (a.nrow (), 1) ;

//			IMatEdit<TTC> resedit = res.Edit () ;

			res.Reset () ;

			DWORD dwAncolm1 = a.ncol () - 1 ;

			DWORD i, h ;
			for (i = a.nrow () - 1; i != (DWORD) -1; i--)
					for (h = dwAncolm1; h != (DWORD) -1; h--)
						res (i, 0) += (TTC) (a (i, h) * b(h)) ;

			return TRUE ;
		}
	} ;


/*
	template <class PAR1, class PAR2>
		class FC_2_0_TYPE
	{
//		friend ;
		protected:
		public:

		static void iter_row_v (const PAR1 &mat, PAR2 &vec, void (*func) (const IVec< typename PAR1::m_datatype> &, typename PAR2::m_datatype &))
		{
			vec.Create (mat.nrow ()) ;
			DWORD i ;
			typename PAR2::m_edittype vecedit (vec.Edit ()) ;
			for (i = mat.nrow (); i != (DWORD) -1; i--)
				func (mat.GetRow (i), vecedit(i)) ;
		}

		static void iter_col_v (const PAR1 &mat, PAR2 &vec, void (*func) (const IVec< typename PAR1::m_datatype > &, typename PAR2::m_datatype &))
		{
			vec.Create (mat.ncol ()) ;
			DWORD i ;
			typename PAR2::m_edittype vecedit (vec.Edit ()) ;
			for (i = mat.ncol (); i != (DWORD) -1; i--)
				func (mat.GetCol (i), vecedit(i)) ;
		}
	} ;

	template <class PAR1, class PAR2>
		class FC_1_1_TYPE
	{
		public:
		static void mean (const PAR1 &vec, PAR2 &val)
		{
			DWORD v ;
			val = 0 ;

			for (v = vec.size () - 1; v != (DWORD) -1; v--)
				val += vec(v) ;
			val /= (PAR2) vec.size () ;
		}
	} ;

	template <class T> inline void mean (const IVec <T> &vec, double &val) { FC_1_1_TYPE<IVec <T>, double>::mean (vec, val) ; }
	template <class T> inline void mean (const IVec <T> &vec, double &val) { FC_1_1_TYPE<IVec <T>, double>::mean (vec, val) ; }

	template <class T> inline double mean (const IVec <T> &vec)	{ double dVal; mean (vec, dVal) ; return dVal ; }
	template <class T> inline double mean (const IVec <T> &vec)	{ double dVal; mean (vec, dVal) ; return dVal ; }
*/
	template <class U, class PAR1, class PAR2>
	void FillMat_Vec_byRow (PAR1 &mat, PAR2 &vec)
	{
//		typename PAR1::m_edittype  matedit = mat.Edit () ;

		const DWORD dwNR = mat.nrow (), dwNC = mat.ncol () ;
		DWORD v = CheckIdx (dwNR * dwNC, vec.size ()) ;

		DWORD c, r ;
		for (c = dwNC - 1; c != (DWORD) -1; c--)
			for (r = dwNR - 1; r != (DWORD) -1; r--)
			{
				mat (c, r) = (typename PAR1::m_datatype) vec (v--) ;
				if (v == (DWORD) -1) v = vec.size () - 1 ;
			}
//		matedit.StopEdit () ;
	}
		
	template <class U, class PAR1, class PAR2>
	void FillMat_Vec_byCol (PAR1 &mat, PAR2 &vec)
	{
//		typename PAR1::m_edittype  matedit = mat.Edit () ;

		const DWORD dwNR = mat.nrow (), dwNC = mat.ncol () ;
		DWORD v = CheckIdx (dwNR * dwNC - 1, vec.size ()) ;

		DWORD r, c ;
		for (r = dwNR - 1; r != (DWORD) -1; r--)
			for (c = dwNC - 1; c != (DWORD) -1; c--)
			{
				mat (c, r) = (typename PAR1::m_datatype)  vec (v--) ;
				if (v == (DWORD) -1) v = vec.size () - 1 ;
			}
	}

#define IMPL_OPERATOR_MAT1(OP, OPAS, FUNC, MTYP, TTYP)																					\
		template <class U> inline IMat<TTYP>  operator OP  (U val)			{ return FC_ElOp<FUNC, TTYP>::OpME	(*this, val); }			\
		template <class U> inline IMat<TTYP> &operator OPAS (U val)			{ return FC_ElOpAs<FUNC>	::OpME	(*this, val); }


/*#define IMPL_OPERATOR_MAT(OP, OPAS, FUNC, MTYP, TTYP)																				\
			template <class U> inline		IMat<TTYP>  operator OP   (const IMat<U> &mat)	const	{ return FC_ElOp<FUNC, TTYP>::OpMM	(*this, mat); }	\
			template <class U> inline const MTYP<TTYP> &operator OPAS (const IMat<U> &mat)	const	{ return FC_ElOpAs<FUNC>	::OpMM	(*this, mat); }	\
			template <class U> inline		IMat<TTYP>  operator OP   (const IVec<U> &vec)	const	{ return FC_ElOp<FUNC, TTYP>::OpMV	(*this, vec); }	\
			template <class U> inline const	MTYP<TTYP> &operator OPAS (const IVec<U> &vec)	const	{ return FC_ElOpAs<FUNC>	::OpMV	(*this, vec); } \
			template <class U> inline		IMat<TTYP>  operator OP   (const U &val)			const	{ return FC_ElOp<FUNC, TTYP>::OpME	(*this, val); }	\
			template <class U> inline const MTYP<TTYP> &operator OPAS (const U &val)			const	{ return FC_ElOpAs<FUNC>	::OpME	(*this, val); }
*/

#define IMPL_OPERATOR_MAT(OP, OPAS, FUNC, MTYP, TTYP)																						\
			inline		IMat<TTYP>  operator OP   (const IMat<TTYP> &mat)		const	{ return FC_ElOp<FUNC, TTYP>::OpMM	(*this, mat); }	\
			inline const MTYP<TTYP> &operator OPAS (const IMat<TTYP> &mat)		const	{ return FC_ElOpAs<FUNC>	::OpMM	(*this, mat); }	\
			inline		IMat<TTYP>  operator OP   (const IVec<TTYP> &vec)		const	{ return FC_ElOp<FUNC, TTYP>::OpMV	(*this, vec); }	\
			inline const	MTYP<TTYP> &operator OPAS (const IVec<TTYP> &vec)	const	{ return FC_ElOpAs<FUNC>	::OpMV	(*this, vec); } \
			inline		IMat<TTYP>  operator OP   (const TTYP &val)				const	{ return FC_ElOp<FUNC, TTYP>::OpME	(*this, val); }	\
			inline const MTYP<TTYP> &operator OPAS (const TTYP &val)			const	{ return FC_ElOpAs<FUNC>	::OpME	(*this, val); }


#define IMPL_OPERATOR_MAT_NOAS(OP, FUNC, MTYP, TTYP)																				\
			template <class U> inline IMat<TTYP>  operator OP   (const U val)			const	{ return FC_ElOp<FUNC, TTYP>::OpME	(*this, val); }	\
			template <class U> inline IMat<TTYP>  operator OP   (const IMat<U> &mat)	const	{ return FC_ElOp<FUNC, TTYP>::OpMM	(*this, mat); }	\
			template <class U> inline IMat<TTYP>  operator OP   (const IVec<U> &vec)	const	{ return FC_ElOp<FUNC, TTYP>::OpMV	(*this, vec); }


#define IMPL_OPERATOR_MAT_BY_ROW(OP, OPAS, FUNC, MTYP, TTYP)																			\
			template <class U> inline IMat<TTYP>  operator OP   (const IVec<U> &vec)	const	{ return FC_ElOp<FUNC, TTYP>::OpMV_row	(Item (), vec); }	\
			template <class U> inline MTYP<TTYP> &operator OPAS (const IVec<U> &vec)		 	{ return FC_ElOpAs<FUNC>	::OpMV_row	(Item (), vec); }
/*
#define IMPL_OPERATOR_VEC(OP, OPAS, FUNC, VTYP, TTYP)																					\
			template <class U> inline const IVec <TTYP> operator OP  (const U val) const	{ return FC_ElOp<FUNC, TTYP>::OpVE	(*this, val); }		\
			template <class U> inline const VTYP<TTYP> &operator OPAS (const U val) const	{ return FC_ElOpAs<FUNC>	::OpVE	(*this, val); }		\
			template <class U> inline const IVec <TTYP>  operator OP (IMat<U> &mat) const	{ return FC_ElOp<FUNC, TTYP>::OpVM	(*this, mat); }		\
			template <class U> inline const IVec <TTYP>  operator OP (IVec<U> &vec) const	{ return FC_ElOp<FUNC, TTYP>::OpVV	(*this, vec); }		\
			template <class U> inline const VTYP<TTYP> &operator OPAS (IVec<U> &vec) const	{ return FC_ElOpAs<FUNC>	::OpVV	(*this, vec); }
*/

#define IMPL_OPERATOR_VEC(OP, OPAS, FUNC, VTYP, TTYP)																					\
			 inline			IVec <TTYP> operator OP  (const TTYP val) const	{ return FC_ElOp<FUNC, TTYP>::OpVE	(*this, val); }		\
			 inline const	IVec<TTYP> &operator OPAS (const TTYP val) const	{ return FC_ElOpAs<FUNC>	::OpVE	(*this, val); }		\
			 inline			IVec <TTYP>  operator OP (IMat<TTYP> &mat) const	{ return FC_ElOp<FUNC, TTYP>::OpVM	(*this, mat); }		\
			 inline			IVec <TTYP>  operator OP (IVec<TTYP> &vec) const	{ return FC_ElOp<FUNC, TTYP>::OpVV	(*this, vec); }		\
			 inline const	IVec<TTYP> &operator OPAS (IVec<TTYP> &vec) const	{ return FC_ElOpAs<FUNC>	::OpVV	(*this, vec); }


#define IMPL_OPERATOR_VEC_NOAS(OP, FUNC, VTYP, TTYP)																					\
			template <class U> inline IVec <TTYP> operator OP  (const U val)	{ return FC_ElOp<FUNC, TTYP>::OpVE	(*this, val); }		\
			template <class U> inline IVec <TTYP>  operator OP (IMat<U> &mat)	{ return FC_ElOp<FUNC, TTYP>::OpVM	(*this, mat); }		\
			template <class U> inline IVec <TTYP>  operator OP (IVec<U> &vec)	{ return FC_ElOp<FUNC, TTYP>::OpVV	(*this, vec); }


