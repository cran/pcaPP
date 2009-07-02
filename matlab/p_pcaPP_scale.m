function x = p_pcaPP_scale ( x )

	x.scl = p_ScaleAdv (x.x, x.center, x.scale) ;
	x.x = x.scl.x ;

	if (isfield (x, 'pHD') && x.pHD) % center and scale must have original data dimension:
		x.scl.center = x.scl.center * (x.svdx.u')  ;
        sca = p_ScaleAdv (x.x * x.svdx.u', null (1), x.scale) ;
		x.scl.scale = sca.scale ;
    end
