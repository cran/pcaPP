#include "pcaPP.h"
#include "qnn.h"

	double ApplyMethod (const SCVecD &v, const int nMethod)
	{
		ASSERT_TEMPRANGE (10, 10) ;
		SVecD temp (tempRef (10), v.size ()) ;	//	2do: should be copied by the constructor!
		temp.Copy_NC (v) ;

		return ApplyMethod_V (!temp, nMethod) ;
	}

	double ApplyMethod_V (const SVVecD &v, const int nMethod)
	{
		double dRet ;
		int nSize = v.size () ;
		switch (nMethod)
		{
			case 0:	sd (dRet, v) ; break ;
			case 1: dRet = mad_V (*v) ; break ;
			case 2: qn (dRet, v.GetData (), nSize) ; break ;
			case 3: dRet = medianabs_V (*v) * 1.482602218505602 ; break ;
//			case 4: sd_st (dRet, v) ; break ;
		}
		return dRet ;
	}
