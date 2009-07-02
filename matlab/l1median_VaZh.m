function [medOut, code, niter] = l1median_VaZh (X, maxiter, tol, zerotol, trace, medIn)

    code = 0 ;
    niter = 0 ;

    if (nargin > 6)
        error ('Too many input arguments.') ;
    elseif (nargin < 6)
        medIn = median (X) ;
        if (nargin < 5)
            trace = 0 ;
            if (nargin < 4)
                zerotol = 1e-15 ;
                if (nargin < 3)
                    tol = 1e-8 ;
                    if (nargin < 2)
                        maxiter = 200 ;
                        if (nargin < 1)
                            error ('Data matrix X is missing.') ;
                        end
                    end
                end
            end
        end
    end

    if (size (X, 1) == 1)
        medOut = X ;
    elseif (size (X, 2) == 1)
        medOut = median (X) ;
    else
        par = [maxiter, trace, tol, zerotol] ;
        dParOut = pcaPP (1, X, medIn, par) ;
        medOut = medIn ;

        code = dParOut(1) ;
        niter = dParOut(2) ;
    end
