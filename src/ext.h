#pragma once

#ifndef _RDEBUG
	#ifdef _MBCS
	#pragma comment(lib, "Rdll.lib") 
	#endif
#endif

#ifdef WIN32
	#define PRE_FUNC extern "C" __declspec(dllexport)
#else
	#define PRE_FUNC extern "C"
#endif //#ifdef WIN32

//#define IMPL_COMP_CONV	//	Compatible Convergence



	PRE_FUNC void l1Median_VZ (unsigned int *pdwPar, double *pdPar, double *pdDat, double *pdMed, double *pdWeights) ;
	PRE_FUNC void l1median_HoCr (int *pnParam_In, int *pnParam_Out, double *pdParam_In, double *pdParam_Out, double *pdData, double *pdMRet) ;

	PRE_FUNC void Hess_Sub_R (int *pnPar, double *pdX_i, double *pdMu, double *pdHess) ;
	PRE_FUNC void Hess_R (int *pnPar, double *pdX, double *pdMu, double *pdHess) ;

	PRE_FUNC void l1median_NM(int *pnParam, double *pdParam, double *pdData, /*double *pdParScale, */double *pdMRet) ;
	PRE_FUNC void l1median_CG(int *pnParam, int *pnParam_Out, double *pdParam, double *pdParam_Out, double *pdData, /*double *pdParScale, */double *pdMRet) ;
	PRE_FUNC void l1median_BFGS (int *pnParam_In, int *pnParam_Out, double *pdParam_In, double *pdParam_Out, double *pdData, /*double *pdParScale, */double *pdMRet) ;
	PRE_FUNC void l1median_SA (int *pnParam_In, int *pnParam_Out, double *pdParam_In, double *pdParam_Out, double *pdData, /*double *pdParScale, */double *pdMRet) ;

	PRE_FUNC void l1median_NLM (int *pnParam, double *pdParam, double *pdData, double *pdMRet/*, double *pdTypSize*/) ;
	PRE_FUNC void l1median_NLM_Hess (int *pnParam, double *pdParam, double *pdData, double *pdMRet/*, double *pdTypSize*/) ;
