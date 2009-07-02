
#ifdef _MSC_VER
	#include "..\..\..\IMat\RFunc.h"
#else
	#include "RFunc.h"
#endif


	double calObj (double *pdData, double *pdM, int n, int p) ;

	class CL1Median_VZ
	{
	public:
		CL1Median_VZ (DWORD *pdwPar, double *pdPar, double *pdDat, double *pdMed, double *pdWeights) ;

		BOOL Iter () ;

		DWORD CheckRowSums (const double &dThreshold) ;

	protected:

		DWORD &m_dwN, m_dwNHalf, &m_dwP, &m_dwMaxIt, &m_dwUseWeights, &m_dwTrace, m_dwEqs, &m_dwTime ;
		double &m_dTol, &m_dZeroTol ;

		IMatD m_mX, m_mX_ ;
		IVecD m_vMed, m_vRt, m_vTt, m_vOldMed ;
		IVecD m_vWeights, m_vRowSums, m_vTemp ;
		IVecB m_mIsZero ;

	} ;
