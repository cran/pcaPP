
#ifdef WIN32
	#define PRE_FUNC extern "C" __declspec(dllexport)
#else
	#define PRE_FUNC extern "C"
#endif //#ifdef WIN32

	PRE_FUNC void rpcgrid (double *px, int *pnParams, double *pObj, double *pLoadings, double *pScores) ;

	PRE_FUNC bool scale_tau (double *pdData, int *pnSize, double dCenter = 0.0, double *pdWeights = NULL, double dInit_Scale = 0.0, double dTuning = 1.95, int bNa_rm = 0) ;
	PRE_FUNC void l1median(double *pdData, int *pnRow, int *pnCol, double *pdMRet, int *pnRet, int *pnMaxStep, double *pdItTol) ;
	PRE_FUNC void l1median_NM(int *pnParam, double *pdParam, double *pdData, double *pdParScale, double *pdMRet) ;
	PRE_FUNC void l1median_CG(int *pnParam, int *pnParam_Out, double *pdParam, double *pdParam_Out, double *pdData, double *pdParScale, double *pdMRet) ;
	PRE_FUNC void l1median_BFGS (int *pnParam_In, int *pnParam_Out, double *pdParam_In, double *pdParam_Out, double *pdData, double *pdParScale, double *pdMRet) ;
	PRE_FUNC void l1median_SA (int *pnParam_In, int *pnParam_Out, double *pdParam_In, double *pdParam_Out, double *pdData, double *pdParScale, double *pdMRet) ;

	PRE_FUNC void l1median_NLM (int *pnParam, double *pdParam, double *pdData, double *pdMRet, double *pdTypSize) ;

	PRE_FUNC void obj (double *pdData, int *pnRows, int *pnCols, double *pdm, double *pds) ;

	PRE_FUNC void rpcnup (double *pdData, int *pnParams, double *pdScores, double *pdVeig, double *pdLambda) ;
	PRE_FUNC void rpcn (double *pdData, int *pnParams, double *pdScores, double *pdVeig, double *pdLambda) ;

	PRE_FUNC void Qn (double *pdX, int *pn, double *pdRet) ;
