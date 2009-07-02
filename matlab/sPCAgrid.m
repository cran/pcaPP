function ret = sPCAgrid (x, k, method, lambda, norm_q, norm_s, maxiter, splitcircle, scores, zero_tol, center, scale, trace)

    if (nargin < 13)
        trace = 0 ;
        if (nargin < 12)
            scale = null (1) ; 
            if (nargin < 11)
                center = 'l1median_HoCr' ;
                if (nargin < 10)
                    zero_tol = 1e-16 ;
                    if (nargin < 9)
                        scores = 1 ;
                        if (nargin < 8)
                            splitcircle = 25 ;
                            if( nargin < 7)
                                maxiter = 10 ;


                                if (nargin < 6)
                                    norm_q = 1 ;
                                    if (nargin < 5)
                                        norm_s = 1 ;

		                        if (nargin < 4)
                                            lambda = 1 ;
                                            if (nargin < 3)
                                                method = 'mad' ;
                                                if (nargin < 2)
                                                    k = 2 ;
                                                    if( nargin < 1)
                                                        error ('Not enough input arguments.') ;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    dat.x = x ;
    dat.k = k ;
    dat.method = method ;
    dat.lambda = lambda ;
    dat.maxiter = maxiter ;
    dat.splitcircle = splitcircle ;
    dat.scores = scores ;
    dat.zero_tol = zero_tol ;
    dat.center = center ;
    dat.scale = scale ;
    dat.trace = trace ;
    dat.ord_all = 0 ;
    

	dat.check_orth = 0 ;

    dat.HDProj = null (1) ;
	dat.HDred = 0 ;
	dat.cut_pc = 1 ;
	dat.glo_scatter = 0 ;
	dat.SpeedUp = 0 ;

	dat = p_sPCAgrid_DataPreProc (dat) ;

    ret_C.nParIn = [dat.k, dat.splitcircle, dat.maxiter, dat.method, dat.trace, dat.k_ini, dat.check_orth, dat.glo_scatter, dat.pHD, dat.SpeedUp] ;
    ret_C.dParIn = [dat.zero_tol, norm_q, norm_s] ;
    ret_C.l = dat.l ;

    [ret_C.nParOut, ret_C.sdev, ret_C.obj] = pcaPP (4 ...
            , ret_C.nParIn ...
            , ret_C.dParIn ...
            , dat.x ...
            , ret_C.l ...
            , dat.lambda ...
            , dat.HDProj ...
    ) ; 

    ret = p_sPCAgrid_PostProc (dat, ret_C) ;
