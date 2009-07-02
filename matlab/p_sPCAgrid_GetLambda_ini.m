function [ lam ] = p_sPCAgrid_GetLambda_ini ( x )

	if (~x.k_ini)
		lam = ones(0) ;
    elseif (is_empty (x.pc_ini.lambda))
        lam = zeros (x.k_ini, 1) ;
    else
        lam = x.pc_ini.lambda ;
    end
