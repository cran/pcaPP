function ret = p_orderPCs ( ret, k, k_ini, ord_all )

    if (nargin < 4)
       ord_all = 0 ; 
    end

	if (isempty (ord_all) || ord_all)
		idx_ord = 1:k ;
	else
		idx_ord = (k_ini + 1):k ;
    end

	if (length (idx_ord) == 1)	%%	nothing to sort_
		return ;
    end

	[soo, ord] = sort (ret.sdev(idx_ord), 'descend') ;

	ord = idx_ord(ord) ;

	ret.pc_order = 1:k ;
	ret.pc_order(idx_ord) = ord ;

	ret.sdev(idx_ord) = ret.sdev(ord) ;
	ret.loadings(:, idx_ord) = ret.loadings (:, ord) ;
	ret.obj(idx_ord) = ret.obj (ord) ;

	if (isfield (ret, 'lambda') && ~isempty (ret.lambda))
		ret.lambda(idx_ord) = ret.lambda (ord) ;
    end
    