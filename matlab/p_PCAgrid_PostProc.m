function [ ret ] = p_PCAgrid_PostProc ( x, ret_C )

	ret.sdev = ret_C.sdev ;
	ret.loadings = ret_C.l ; %matrix (ret_C.l, ncol = x.p, nrow = x.p)
    ret.k = x.k ;
    ret.obj = ret_C.obj ;
    ret.n_obs = size (x.x_orig, 1) ;
    %ret.args = x.args ;
	ret.scale = x.scl.scale ;
    ret.center = x.scl.center ;

	if(x.pHD)			%%	undo SVD for high dimensional datasets
		ret.loadings = x.HDProj * ret.loadings ;
    end

	if (x.cut_pc)
		ret = p_cut_pc (ret, x.k) ;
    end

	ret.loadings = p_loadSgnU (ret.loadings) ;

	ret = p_orderPCs (ret, x.k, x.k_ini, x.ord_all) ;

    
    if (x.scores)
        ox = ones (size (x.x_orig, 1), 1) ;
		ret.scores = (x.x_orig - ox * x.scl.center) * ret.loadings ;
    end
