function ret = p_cut_pc (ret, k)

	ret.loadings = ret.loadings (:, 1:k) ;
	ret.sdev = ret.sdev (1:k) ;
	if (isfield (ret, 'scores') && ~isempty (ret.scores))
		ret.scores = ret.scores(:, 1:k) ;
    end
	if (isfield (ret, 'lambda') && ~isempty (ret.lambda))
		ret.lambda = ret.lambda (1:k) ;
    end
