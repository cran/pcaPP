function [ ret ] = p_sPCAgrid_PostProc (x, ret_C)

	ret = p_PCAgrid_PostProc (x, ret_C) ;
    
    lambda_ini = p_sPCAgrid_GetLambda_ini (x) ;

	ret.lambda = [lambda_ini, x.lambda] ;
