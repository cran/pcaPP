function [ ret ] = qn ( x, corrfact )

    if (nargin < 2)
        corrfact = 2.219144465985075420633 ; % 1/(sqrt(2) * qnorm(5/8))
        if (nargin < 1)
             error ('Not enough input arguments.')
        end
    end

    ret = pcaPP (2, x, corrfact) ;
