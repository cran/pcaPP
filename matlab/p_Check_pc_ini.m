function x = p_Check_pc_ini (x)
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here

	if (isfield (x,'pc_ini'))
		if (~isfield (x, 'k_ini'))
			x.k_ini = x.pc_ini.k ;
        end
	else
		x.k_ini = 0 ;
    end        

	if (x.k_ini)

        if (any (size (x.pc_ini.load) ~= size (x.x, 2)))
            error ('pc_ini.load has wrong dimension.')
        end
        if (x.k_ini + x.k > size (x.x, 2))
            error ('k has been chosen too large') ;
        end
		if (isempty (x.k_ini))
   			x.k_ini = x.pc_ini.k ;
        end
        
		x.k = x.k_ini + x.k	;%%	2do: don't change k anymore to k_ini + k
        if (length (x.pc_ini.sdev) ~= ncol (x.x))
            error ('vector pc_ini.sdev has wrong size') ;
        end
		
		x.sdev = x.pc_ini.sdev ;
		x.l = x.pc_ini.load ;
	else
		x.k_ini = 0 ;
        x.sdev = ones (x.p, 1) * NaN ;
        x.l = diag (ones (x.p, 1)) ;
    end



