#include "PCAgrid.h"

#define CHECK_ORTH
#define GLOBAL_POWER	2	//	use 2 for maximizing (robust) variances, 1 for maximizing (robust) standard deviations.

#if GLOBAL_POWER == 1
	const &double ngpf (const double &d) {return d ; }
#else 
	#if GLOBAL_POWER == 2
		double ngpf (const double &d) {return d * d ; }
	//	#else: error, because gpf not defined.
	#endif
#endif

/////////////////
//	CsPCAGrid  //
/////////////////


	CsPCAGrid::CsPCAGrid (int *pnParamIn, int *pnParamOut, double *pdParamIn, double *pdData, double *pdLoadings, double *pdSDev, double *pdObj/*, double *pdMaxMaha*/, double *pdLambda, double *pdBackTransHD)
		:	CPCAGrid (pnParamIn, pnParamOut, pdParamIn, pdData, pdLoadings, pdSDev, pdObj/*, pdMaxMaha*/)
		,	m_nGloScatter (pnParamIn[9]), m_nSpeedUp (pnParamIn[11]), m_dwPHD (pnParamIn[10])
		,	m_vLambda (pdLambda, m_dwK), m_vTempP (m_dwP), m_vTempPSub (m_dwP)
		,	m_dGloScatter (1)
	{
		if (m_dwPHD)
		{
			m_mBackTransHD.Set (pdBackTransHD, m_dwPHD, m_dwP) ;
			m_mBackProj.Require (m_dwP, m_dwPHD) ;
			m_vLoadHD.Require (m_dwPHD) ;
			m_vSumLoadOthers.Require (m_dwPHD) ;
		}
		else
		{
			m_mBackProj.Require (m_dwP, m_dwP) ;
			m_vSumLoadOthers.Require (m_dwP) ;
		}
		if (!m_dwCheckOrth && m_nGloScatter == 0)
			m_dGloScatter = ngpf (ApplyMethodMean (m_mX)) ;
//		else if (m_nGloScatter == -1)
//			m_dGloScatter = 1 ;
	}

	void CsPCAGrid::OnCalcPC ()
	{
		if (!m_dwCheckOrth && m_nGloScatter == 1)
			m_dGloScatter = ngpf (ApplyMethodMean (TempY ())) ;
//meal_printf ("gloscat %f\r\n", m_dGloScatter) ;

		m_vTempPSub.Reshape (m_dwPSub) ;
		m_dCurLambda = m_vLambda (m_dwCurK - m_dwkIni) ;

//meal_printf ("Glo scat for dim %d = %f\r\n", m_dwCurK, sqrt (m_dGloScatter)) ;

		// calc new backtransformation matrix.

		if (m_dwPHD)
			sme_matmult_R (m_mBackTransHD, m_mL.GetColRef (m_dwCurK, m_dwP), !m_mBackProj) ;
		else
			m_mBackProj.Copy_R (m_mL.GetColRef (m_dwCurK, m_dwP)) ;	//	2do: use a reference!
	}

	void CsPCAGrid::InitPenalty ()
	{
		m_vSumLoadOthers.Reset (0) ;
		EO<SOP::ApaBmC>::VMcVct (*m_vSumLoadOthers, m_mBackProj, m_vAfin) ;

		m_vSumLoadThis.Copy_R (m_mBackProj.GetColRef (m_dwCurP)) ;	//	2do: try simple assignment or use this column directly where m_vSumLoadThis is needed!
	}

	double CsPCAGrid::CalcObj (const double dCos, const double dSin, double &dScat, double &dScatOrth)
	{
		dScat = CalcProjScat (dCos, dSin) ;

		if (m_dwCheckOrth)
		{
//			CalcMaha (dScat) ;


			dScatOrth = CalcProjScat (-dSin, dCos) ;
//			CalcMaha (dScatOrth) ;

			return ngpf (dScat / dScatOrth) + GetPenalty (dCos, dSin); ;

////			dCurObj /= ngpf (dScatOrth + dCurScat * 0.001) ;
////			dCurObj /= (m_dGloScatter + ngpf (dScatOrth)) ;
		}
		return ngpf (dScat) + GetPenalty (dCos, dSin) * m_dGloScatter;
	}

	double CsPCAGrid::GetPenalty (const double& dCos, const double& dSin)
	{
		ASSERT (!m_nSpeedUp) ;	//	the old penalty stuff has not been implemented here...

		if (m_dCurLambda == 0)
			return 0 ;

		double dRet = 0 ;

		if (fabs (dCos) <= m_dZeroTol)
			EO<UOP::Apa_abs_B>::SVc (dRet, m_vSumLoadThis) ;
		else if (fabs (dSin) <= m_dZeroTol)
			EO<UOP::Apa_abs_B>::SVc (dRet, m_vSumLoadOthers) ;
		else
			EO<UOP::Apa_abs_BmDpCmE_>::SScScVcVc_NC (dRet, dCos, dSin, m_vSumLoadOthers, m_vSumLoadThis) ;

		return - dRet * m_dCurLambda ;

//		dRet *= - m_dCurLambda ;
//		if (!m_dwCheckOrth)
//			dRet *= m_dGloScatter ;
//
//		return  dRet ;
	}

////////////////
//	CPCAGrid  //
////////////////
						

	CPCAGrid::CPCAGrid (int *pnParamIn, int *pnParamOut, double *pdParamIn, double *pdData, double *pdLoadings, double *pdSDev, double *pdObj/*, double *pdMaxMaha*/)
		: m_dwN (pnParamIn[0]), m_dwP (pnParamIn[1]), m_dwK (pnParamIn[2]), m_dwSplitCircle (pnParamIn[3]), m_dwMaxIter (pnParamIn[4]), m_dwMethod (pnParamIn[5]), m_dwTrace (pnParamIn[6]), m_dwkIni (pnParamIn[7]), m_dwCheckOrth (pnParamIn[8])
		, m_nReturn (pnParamOut [0])
		, m_dZeroTol (pdParamIn[0])
		, m_mX (pdData, m_dwN, m_dwP), m_mL (pdLoadings, m_dwP, m_dwP)//, m_mTempPP(m_dwP, m_dwP), m_mTempPN(m_dwN, m_dwP)
		, m_vAfin(m_dwP), m_vAfinBest (m_dwP), m_vScl(m_dwP), m_vYOpt (m_dwN), m_vSDev(pdSDev, m_dwP), m_vObj (pdObj, m_dwK)
		, m_vProj (m_dwN)//, m_vMaxMaha (pdMaxMaha, m_dwN)
		, m_vOrd (m_dwP)
		, m_dwCurK (0), m_dwCurP (0), m_dwPSub (0), m_dwTempYIdx (0)

		, m_pdProj (m_vProj)//, m_pdEndProj (m_vProj.GetDataEnd ())
			, m_pdCurLC (m_vYOpt), m_pdCurLCEnd (m_vYOpt.GetDataEnd ())
	{
		m_mY[0].Require (m_dwN, m_dwP) ;
		m_mY[1].Require (m_dwN, m_dwP) ;
//		m_vMaxMaha.Reset (0) ;
	}

	int CPCAGrid::Calc ()
	{
		if (m_dwK > m_dwP)
			return 1 ;		//	k > p

/*		if ((m_nSplitCircle & 1) == 0)	//
		//if (m_dwSplitCircle & 1)		//  only allow even values for splitcircle
			++ m_nSplitCircle ;*/

		if (m_dwkIni)
			sme_matmult_R (m_mX, m_mL.GetColRef (m_dwkIni, m_dwP), !TempY ()) ;
		else
		{
			TempY ().Copy (m_mX) ;
			SetDiag_sq (!m_mL) ;
			//m_mL.setdiag () ;			//	this MUST now happen in the calling R routine! // this has been changed when introducing the m_dwkIni argument // why not here?
		}

		for (m_dwCurK = m_dwkIni; m_dwCurK < m_dwK; m_dwCurK++)
		{
			m_dwPSub = m_dwP - m_dwCurK ;		//	dimensionality of the subspace

			OnCalcPC () ;

			if (m_dwPSub == 1)
			{
				m_vSDev (m_dwCurK) = ApplyMethod (TempY ().GetColRef (0)) ;
				continue ;	//	break ;
			}

			m_vScl.Reshape (m_dwPSub) ;
			m_vOrd.Reshape (m_dwPSub) ;
			ApplyMethod (TempY (), m_vScl) ;	//	2do: m_vScl can be a temporary vector

			meal_sort_order_rev (m_vScl, m_vOrd, m_vScl.size ()) ;

			m_dwCurP = m_vOrd(0) ;				//	index of the coloumn of x with biggest scatter

			m_vAfinBest.Reshape (m_dwPSub) ;
			m_vAfin.Reshape (m_dwPSub) ;

			m_vAfin.Reset (0) ;
			m_vAfin (m_dwCurP) = 1 ;

			CopyCol (*m_vYOpt, TempY (), m_dwCurP) ;

			t_size i, j ;
			double dCurSplit ;

			double dScatBest = 0 ;

			double dObjBest = 0 ;
			for (i = 0; i <= m_dwMaxIter; i++)
			{
				//double dScat, dObj, dSumAbsDelta = 0 ;

				dCurSplit = pow (0.5, (double) i) ;		// or use  * 0.5 each time?

				for (j = 0; j < m_dwPSub; j++)
				{
					m_dwCurP = m_vOrd (j) ;
					m_vCurY = TempY ().GetColRef (m_dwCurP) ;	//	2do: move this 2 rows down.
					m_pdCurY = m_vCurY ;

					const double dL = m_vAfin (m_dwCurP) ;	//	current loading 
					if (fabs (dL) == 1)
						continue ;
					RemoveLoading (/*i*/) ;

					m_dNL = dL ;
					GridPlane (dCurSplit) ;
					AddLoading (m_dNL, m_dNCL) ;

					//double dNL = dL, dNCL ;
					//GridPlane (dNL, dNCL, dScat, dObj, dCurSplit) ;
					//AddLoading (dNL, dNCL) ;
					//dSumAbsDelta += fabs (dL - dNL) ;
				}

				EO<SOP::a_divide>::VSc (*m_vAfin, norm2 (m_vAfin)) ;	//	2do: check norm of m_vAfin. should be 1 anyway!, if not it's sufficient to perform this normalization after this for loop, only with the m_vAfinBest - vector!

				if (!i || dObjBest <= m_dBestObj)
				{
					dObjBest = m_dBestObj ;
					m_vAfinBest.Copy_NC (m_vAfin) ;
					dScatBest = m_dCurScat ;
				}

//meal_printf ("delta: %.22f ->", dSumAbsDelta) ;
/*				if (dSumAbsDelta <= m_dZeroTol)	//	no changes of any loading -> quit // doesn't make sense as we're operating on a raster
				//if (dCurSplit<= m_dZeroTol)	//	no changes of any loading -> quit
				{
//meal_printf ("stop iteration\n") ;
					if (m_dwTrace >= 3)
						meal_printf ("Calculation of PC %d stopped after %d loops\r\n", m_dwCurK + 1, i + 1) ;
					break ;
				}
//meal_printf ("continue iteration\n") ;
*/			}

			m_vSDev (m_dwCurK) = dScatBest ;		//	2do: use ptrs instead of vector access!
			m_vObj(m_dwCurK) = dObjBest ;
			BackTransform () ;
		}
		return 0 ;
	}

	void CPCAGrid::BackTransform ()
	{
		ASSERT_TEMPRANGE (0, 2) ;
		SMatD m_mTempPP (tempRef (0), m_dwPSub, m_dwPSub) ;
		SetDiag_sq (!m_mTempPP) ;

		int dwIdxRef = m_vOrd (0) ;

		set_neg (*m_vAfinBest) ;

		m_vAfinBest (dwIdxRef) += 1 ;

		double dNorm = norm2 (m_vAfinBest) ;

		if (dNorm > m_dZeroTol)
		{
			static const double dSqrt2 = sqrt ((double) 2.0) ;
			EO<SOP::a_divide>::VSc (*m_vAfinBest, dNorm / dSqrt2) ;
			EO<SOP::AsaBmC>::MVcVct (!m_mTempPP, m_vAfinBest, m_vAfinBest) ;
		}

		SMatD mProjSorted (tempRef (1), m_dwPSub, m_dwPSub) ;

		mProjSorted.CopyCol_Order_NC (m_mTempPP, *m_vOrd) ;						//	undo some kind of sorting

		SMatD mOldLoadings (tempRef (2), m_dwP, m_dwPSub) ;						//	2do: copying cols should be done in constructor (using GetColRef ()
		CopyCol (!mOldLoadings, m_mL, m_dwCurK, m_dwP) ;

		sme_matmult (mOldLoadings, mProjSorted, !m_mL.GetColRef (m_dwCurK, m_dwP)) ;
		sme_matmult_R (TempY (), mProjSorted.GetColRef (1, m_dwPSub), !TempYC ()) ;

		SwapTempY () ;
	}

	void CPCAGrid::RemoveLoading (/*int i*/)
	{
		const double &dL = m_vAfin (m_dwCurP) ;
		if (dL == 0)
			return ;

		const double dCL = sqrt (1.0 - sm_sqr (dL)) ;	//	current loading and complementary loading, such that dL^2 + dCL^2 == 1

		EO<UOP::Aa_AsDmB_dC>::VScScVc (*m_vYOpt, dL, dCL, m_vCurY) ;
		EO<SOP::a_divide>::VSc (*m_vAfin, dCL) ;
		m_vAfin(m_dwCurP) = 0 ;
	}

	void CPCAGrid::AddLoading (const double &dNL, const double &dNCL)
	{
		//	dNL		...	New Loading
		//	dNCL	... New Complement Loading

		EO<UOP::Aa_AmC_p_DmB>::VScScVc (*m_vYOpt, dNL, dNCL, m_vCurY) ;

		EO<SOP::a_multiply>::VSc (*m_vAfin, dNCL) ;
		m_vAfin (m_dwCurP) = dNL ;
	}

	double CPCAGrid::ApplyMethodMean (const SCMatD &m)
	{
		double dSd = 0 ;
		int i ;
		for (i = m.ncol () - 1; i != (int) -1; i--)
			dSd += sm_sqr(ApplyMethod (m.GetColRef (i))) ;
		return sqrt (dSd  / m.ncol ()) ;
	}

	void CPCAGrid::ApplyMethod (const SCMatD &m, SVecD &v)
	{
		v.Reshape (m.ncol ()) ;
		int i ;
		for (i = m.ncol () - 1; i != (int) -1; i--)
			v(i) = ApplyMethod (m.GetColRef (i)) ;
	}

	double CPCAGrid::ApplyMethod (const SCVecD &v)	//	2do: remove..
	{
		return ::ApplyMethod (v, m_dwMethod) ;
	}

	double CPCAGrid::CalcProjScat (const double dCos, const double dSin)
	{
		const double *pYOpt = m_pdCurLC, *pCurY = m_pdCurY ;
		double *pProj = m_pdProj ;

		while (pYOpt < m_pdCurLCEnd)
		{
			*pProj = *pYOpt * dCos + *pCurY  * dSin ;
			++pProj ;
			++pYOpt ;
			++pCurY ;
		}

		return ApplyMethod (m_vProj) ;
	}

/*	void CPCAGrid::CalcMaha (const double dScat)
	{
		double *pProj = m_pdProj ;
		double *pdProjEnd = m_pdProj + m_dwN ;

		double *pdMaha = m_vMaxMaha ;

		while (pProj < pdProjEnd)
		{
			sm_setmax (*pdMaha, *pProj / dScat) ;
			++pdMaha ;
			++pProj ;
		}
	}
*/
	double CPCAGrid::CalcObj (const double dCos, const double dSin, double &dScat, double &dScatOrth)
	{
		dScat = CalcProjScat (dCos, dSin) ;

		if (m_dwCheckOrth)
		{
//			CalcMaha (dScat) ;

			dScatOrth = CalcProjScat (dCos, -dSin) ;
//			CalcMaha (dScatOrth) ;

			return ngpf (dScat / dScatOrth) ;

////			dCurObj /= ngpf (dScatOrth + dCurScat * 0.001) ;
////			dCurObj /= (m_dGloScatter + ngpf (dScatOrth)) ;
		}

		return ngpf (dScat) ;
	}

	void CPCAGrid::EvalDirection (const double dCos, const double dSin)
	{
		double	dScat, dScatOrth ;
		double dCurObj = CalcObj (dCos, dSin, dScat, dScatOrth) ;

		if (dCurObj > m_dBestObj)
		{
			m_dBestObj = dCurObj ;
			m_dCurScat = dScat ;
			m_dCurScatOrth = dScatOrth ;

			m_dNL = dSin ;
			m_dNCL = dCos ;
		}
	}

	double CPCAGrid::CalcVarTrimmed (double dCos, double dSin, double dScat, double dScatOrth)
	{
		if (dScatOrth <= m_dZeroTol || dScat <= m_dZeroTol)
			return dScat ;

		const double *pYOpt = m_pdCurLC, *pCurY = m_pdCurY ;
		double dCurProj, dCurProjOrth ;

		dScat = 1 / dScat ;
		dScatOrth = 1 / dScatOrth ;

		double dS = 0, dSS = 0 ;
		int n = 0 ;

		while (pYOpt < m_pdCurLCEnd)
		{
			dCurProj = *pYOpt * dCos + *pCurY  * dSin ;
			dCurProjOrth = *pYOpt * dSin - *pCurY  * dCos ;

			if (sm_sqr (dCurProj) * dScat + sm_sqr (dCurProjOrth) * dScatOrth < 6)
			{
				dS += dCurProj ;
				dSS += sm_sqr (dCurProj) ;
				++n ;
			}

			++pYOpt ;
			++pCurY ;
		}

		return (dSS / n - (sm_sqr (dS / n))) * n / (n - 1.0) * 1.3178;	//	correction factor 1.3178 for 95% quantile in sqr maha distance (6)
	}

	double CPCAGrid::CalcScatTrimmed (double dCos, double dSin, double dScat, double dScatOrth)
	{
		if (dScatOrth <= m_dZeroTol || dScat <= m_dZeroTol)
			return dScat ;

		const double *pYOpt = m_pdCurLC, *pCurY = m_pdCurY ;
		double dCurProjOrth ;
		double *pProj = m_pdProj ;

		while (pYOpt < m_pdCurLCEnd)
		{
			dCurProjOrth = *pYOpt * dSin - *pCurY  * dCos ;

			if (sm_sqr (dCurProjOrth) / dScatOrth <= 3.841459)
			{
				*pProj = *pYOpt * dCos + *pCurY  * dSin ;
				++pProj ;
			}

			++pYOpt ;
			++pCurY ;
		}

//		return (dSS / n - (sm_sqr (dS / n))) * n / (n - 1.0) * 1.3178;	//	correction factor 1.3178 for 95% quantile in sqr maha distance (6)
		return ApplyMethod (SVecD (m_pdProj, pProj - m_pdProj)) ;
	}

	void CPCAGrid::GridPlane (double dCurSplit)
	{
		const double dASinNL = asin (m_dNL) ;

		const double dSm1 = (m_dwSplitCircle > 1) ? m_dwSplitCircle - 1 : 1 ;

		double dSplitFact = meal_PI () * dCurSplit ;

		InitPenalty () ;

		t_size i ;

		ASSERT_TEMPRANGE (11, 11) ;

		SVecD vProj (tempRef (11), m_dwN) ;

		m_dBestObj = meal_NegInf () ;

		if (m_dNL && fabs (m_dNL) < 1e-6)
			EvalDirection (1, 0) ;

		const t_size dwEnd = (dCurSplit == 1.0) ? m_dwSplitCircle - 1 : m_dwSplitCircle ;	//	dCurSplit means, that we're checking an angle of 180° ( = PI). thus the first and last checked point would be the same.

		for (i = 0; i < dwEnd; i++)
		{
			double dAngle = (i / dSm1 - 0.5) * dSplitFact + dASinNL;

			EvalDirection (cos (dAngle), sin (dAngle)) ;
		}
		if (m_dwCheckOrth)
			m_dCurScat = sqrt (CalcVarTrimmed (m_dNCL, m_dNL, m_dCurScat, m_dCurScatOrth)) ;

//			m_dCurScat = CalcScatTrimmed (m_dNCL, m_dNL, m_dCurScat, m_dCurScatOrth) ;
//		else
//			m_dCurScat = m_dCurScat ;
	}
