#include "L1Median.h"
#include "ext.h"
#include "perftimer.h"

	DWORD CL1Median_VZ::CheckRowSums (const double &dThreshold)
	{		//	counts the elements in array m_vRowSums which are greater than dThreshold
		DWORD i ;
		DWORD dwRet = 0 ;
		for (i = m_vRowSums.size () - 1; i != (DWORD) -1; i--)
			if (m_mIsZero(i) = m_vRowSums  (i) > dThreshold)
				dwRet ++ ;
		return m_dwN - dwRet ;	
	}

	BOOL CL1Median_VZ::Iter ()
	{
		m_mX_ << m_mX ;
		m_mX_.byrow () -= m_vMed ;	//	centering data matrix

		SqrtrowSumSqs (m_mX_, m_vRowSums) ;

		double dMin = min (m_vRowSums) ;

		DWORD dwGreater = m_dwN - m_vRowSums.CountLess (dMin / m_dZeroTol) ;

		if (dwGreater * 2 > m_dwN)
		{	//	some of the elements of m_vRowSums are zero
			m_dwEqs++ ;

			//double dCheckMedian = median (m_vRowSums, m_vTemp) ;
			DWORD dwZero = CheckRowSums (median (m_vRowSums, m_vTemp) * m_dZeroTol) ;

			if (dwZero > m_dwNHalf)	//	there's more than half of the observations concentrated in one single point
			{
				if (m_dwTrace >= 1)
					Rprintf ("%d >= n / 2 = %d observations concentrated in one point found.\r\n", dwZero) ;
				return FALSE ;
				
			}

			if (m_dwTrace >= 1)
				Rprintf ("%d observations are exatly at the median.\r\n", dwZero) ;
			if (int (m_dwTrace) >= 0 && dwZero > 1)
				Rf_warning ("The current L1median estimate is ident with more than one observation. The resulting l1median estimation might be incorrect. [CL1Median_VZ::Iter]") ;

			IMatD	X_0	(m_mX_	(ISub (m_mIsZero), ISub ())), 
					X0	(m_mX	(ISub (m_mIsZero), ISub ())) ;
			IVecD	vID (m_vRowSums (ISub (m_mIsZero))) ;

			vID ^= -1 ;

			ColSumWeighted (X_0, vID, m_vRt) ;
			ColSumWeighted (X0, vID, m_vTt) ;
			m_vTt /= sum (vID) ;

			double edivr = dwZero / sqrtsumsq (m_vRt) ;
			if (edivr < 1)
				m_vMed *= edivr ;

			if(edivr < 1)
			{
				m_vTt *= 1 - edivr ;
				m_vMed += m_vTt ;
			}
		}
		else
		{
			m_vRowSums ^= -1 ;

			ColSumWeighted (m_mX, m_vRowSums, m_vMed) ;
			m_vMed /= sum (m_vRowSums) ;
		}
		return TRUE ;
	}

	CL1Median_VZ::CL1Median_VZ (DWORD *pdwPar, double *pdPar, double *pdDat, double *pdMed, double *pdWeights) :
		m_dwN (pdwPar[0]), m_dwNHalf (m_dwN >> 1), m_dwP (pdwPar[1]), m_dwMaxIt (pdwPar[2]), m_dwUseWeights (pdwPar[3]), m_dwTrace (pdwPar[4]), m_dwEqs (0), m_dwTime (pdwPar[5]),
		m_dTol (pdPar[0]), m_dZeroTol (pdPar [1]),
		m_mX (m_dwN, m_dwP, pdDat), m_mX_ (m_dwN, m_dwP),
		m_vMed (m_dwP, pdMed), m_vRt (m_dwP), m_vTt (m_dwP), m_vOldMed (m_dwP), 
		m_vWeights (m_dwN, pdWeights), m_vRowSums (m_dwN), m_vTemp (m_dwN),
		m_mIsZero (m_dwN)
	{ 
		DWORD i ;

		double dAbsDiff, dAbsSum ; 

		CPerfTimer tim ;
		for (i = m_dwMaxIt - 1; i != (DWORD) -1; i--)
		{
			m_vOldMed << m_vMed ;
			if (!Iter ())
				break ;
			m_vOldMed -= m_vMed ;

			dAbsDiff = sumabs (m_vOldMed) ;
			dAbsSum = sumabs (m_vMed) ;

			if (m_dwTrace >= 2)
				if (m_dwTrace >= 3)
				{
					Rprintf ("k=%3d rel.chg=%12.15g, m=(",m_dwMaxIt - i, dAbsDiff/dAbsSum),
					Rprintf (")\n") ;
				}
				else
					Rprintf (".") ;


			if (dAbsDiff < m_dTol * dAbsSum)
				break ;
		}
		m_dwTime = (DWORD) tim.GetTimeMS () ;

        if(m_dwTrace)
			Rprintf (" needed %d iterations (%d of which had y==x_k)\r\n", m_dwMaxIt - i, m_dwEqs) ;

		m_dwMaxIt -= i ;

	}

	
	
	EXPORT void l1Median_VZ (DWORD *pdwPar, double *pdPar, double *pdDat, double *pdMed, double *pdWeights)
	{
		CL1Median_VZ (pdwPar, pdPar, pdDat, pdMed, pdWeights) ;
	}

