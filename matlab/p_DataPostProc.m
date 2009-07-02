function [ret] = p_DataPostProc (DataObj, obj, loadings, scores, bScores)
    
    %idx = order (obj, decreasing = true) ;
    [soo, idx] = sort (obj, 'descend') ;

    obj = obj (idx) ;
    loadings = loadings (:, idx) ;

	if (bScores)
		scores = scores (:, idx) ;
    end

%	ret = list()

   %%loadings
	
		ret.loadings = loadings ;
		ret.loadings = p_loadSgnU (ret.loadings) ;
        
%		c = ncol (loadings)
%		r = nrow (loadings)
%		if (is.null (dimnames (DataObj.x)[[2]]))
%			dimnames (ret.loadings) = list (paste (rep ('V', r), 1:r, sep = ''), paste (rep ('Comp.', c), 1:c, sep = ''))
%		else
%			dimnames (ret.loadings) = list (dimnames (DataObj.x)[[2]], paste (rep ('Comp.', c), 1:c, sep = ''))

%		class (ret.loadings) = 'loadings'

   %%sdev
	ret.sdev = obj ;
	%names (ret.sdev) = dimnames (ret.loadings)[[2]]

   %%center
	ret.center = DataObj.center ;
   %%scale
	ret.scale = DataObj.scale ;
   %%n_obs
	ret.n_obs = size (DataObj.x, 1) ;

   %%scores
	if (bScores)
		ret.scores = scores ;
		%dimnames (ret.scores) = list (1:nrow (scores), dimnames (ret.loadings)[[2]]) ;
	else
		ret.scores = null (1) ;
    end
%	ret.call = cl

%	class (ret) = c ('pcaPP', 'princomp')
%	return (ret)
%}
