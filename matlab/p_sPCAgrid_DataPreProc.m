function x = p_sPCAgrid_DataPreProc ( x )
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here



	len_lambda = length (x.lambda) ;

	if (len_lambda ~= 1 && len_lambda ~= x.k)
		warning ('length (lambda) should either be equal to 1 or k') ;
    end
	x.lambda = p_replen (x.lambda, x.k) ;

	x = p_PCAgrid_DataPreProc (x) ;
