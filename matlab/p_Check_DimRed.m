function x = p_Check_DimRed ( x )

    x.pHD = 0 ;
    
    x.x_orig = x.x ;

    if (x.p > x.n && x.HDred)			%% Dimension reduction for high dimensional datatsets

		[svdx.u, S, svdx.v] = svd(x.x') ;
        
        svdx.u = svdx.u (:, 1:size(x.x,1)) ;
        svdx.d = diag (S) ;
        
		x.svdx = svdx ;

		x.x = svdx.v * diag (svdx.d) ;
		x.HDProj = svdx.u ;

		x.pHD = x.p ;
		[x.n, x.p] =  size (x.x) ;

		if (x.trace >= 2)
            disp (sprintf ('reduced dimensions -> n x p = %dx$x\n', x.n, x.p))
        end
    end

