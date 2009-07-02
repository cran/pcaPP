function [ ret ] = PCAproj (x, k, method, CalcMethod, nmax, update, scores, maxit, maxhalf, scale, center, zero_tol)

    if (nargin < 12)
        zero_tol = 1e-16 ;
        if (nargin < 11)
            center = 'l1median_HoCr' ;
            if (nargin < 10)
                scale = null (1) ;
                if (nargin < 9)
                    maxhalf = 5 ;
                    if (nargin < 8)
                        maxit = 5 ;
                        if (nargin < 7)
                            scores = true ;
                            if (nargin < 6)
                                update = true ;
                                if (nargin < 5)
                                    nmax = 1000 ;
                                    if (nargin < 4)
                                        CalcMethod = 'eachobs' ;
                                        if (nargin < 3)
                                            method = 'mad' ;
                                            if (nargin < 2)
                                                k = 2 ;
                                                if (nargin < 1)
                                                    error ('Not enough input arguments provided.') ;
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
    
	method = p_GetScaleMethod (method) ;

	[n, p] = size (x) ;

	if( k > min(n,p))
		error ('k too large')
    end

	if(p > n)
        
            [svdx.u, S, svdx.v] = svd(x') ;

            svdx.u = svdx.u (:, 1:size(x,1)) ;
            svdx.d = diag (S) ;
    		x = svdx.v * diag (svdx.d) ;

%			svdx = svd(t(x)) ;
%			x = svdx.v * diag(svdx.d) ;

			pold = p ;
			p = n ;
	else
			pold = p ;
    end

	DataObj = p_ScaleAdv (x, center, scale) ;

    if (pold > n) % center and scale must have original data dimension:
		DataObj.center = svdx.u * DataObj.center' ;
        scaled = p_ScaleAdv (x * svdx.u', null (1), scale) ;
		DataObj.scale = scaled.scale ;
    end

	y = DataObj.x ;

    if (scores)
		scoresize = n * k ;
	else
		scoresize = 0 ;
    end

	if (strcmp (CalcMethod, 'lincomb'))
        update = false ;
		if (nmax > n)
            aux = rand (nmax - n, n) ;
			%aux = matrix (runif ((nmax-n) * n), nrow = nmax-n) ;
            yadd = (aux * x) - (ones (size (aux, 1), 1) * DataObj.center) ;
            %yadd = (aux * y) ; % shouldn't this be better?
            y = [y ; yadd] ;
			%y = rbind (y, t(t(aux * as.matrix(x)) - DataObj.center)) ;
        end
    elseif (strcmp (CalcMethod, 'sphere'))
        update = false ;
		if(nmax > n)
            yadd = randn (nmax - n, p) ;
            y = [y ; yadd] ;
			% y = rbind (y, randn (nmax-n, p)) ; %rmvnorm(nmax-n, rep(0,p), diag (p))) ;
        end
    elseif  (strcmp (CalcMethod, 'eachobs'))
    else
       error ('Argument CalcMethod: invalid value.')
    end

	nn = size (y, 1) ;

	if (update)
        nParIn = [n, k, method, scores, maxit, maxhalf] ;
        [z, l, s] = pcaPP (6, nParIn, zero_tol, y) ;
    else
        nParIn = [n, k, method, scores] ;
        [z, l, s] = pcaPP (5, nParIn, zero_tol, y) ;
    end

	%veig = matrix (ret.C.loadings, ncol = k) ;

    if (pold > n)
		l = svdx.u * l ;
    end

	if (scores)
		ret = p_DataPostProc (DataObj, s, l, z(1:n, :) , scores) ;
	else
		ret = p_DataPostProc (DataObj, s, l, null (1), scores) ;
    end
