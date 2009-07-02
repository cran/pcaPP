function x = p_PCAgrid_DataPreProc ( x )
%mydesc
	%%	checking/initializing parameters

	[x.n, x.p] = size (x.x) ;

    if (~isnumeric (x.trace))
        error ('trace must be numeric.') ;
    end
    if (length (x.k) ~= 1)
        error ('k must be a scalar.') ;
    end

    if (x.k < 1)
        error ('k must be larger than 0.') ;
    end

	x.method = p_GetScaleMethod (x.method) ;

	x = p_Check_DimRed (x) ;
	x = p_Check_pc_ini (x) ;

    if (x.k > size (x.x, 2))
       error ('k has been chosen too large.') 
    end

	x = p_pcaPP_scale (x) ;

