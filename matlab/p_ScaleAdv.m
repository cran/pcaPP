function [ ret ] = p_ScaleAdv ( x, center, scale )

    if (nargin < 3)
       scale = null (1) ;
       if (nargin < 2)
          center = null  (1) ;
          if (nargin < 1)
             error ('Not enough input arguments.')
          end
       end
    end

    p = size (x, 2) ;

    ret.center = p_genapply(x, center, 'center', 0) ;
    ret.scale = p_genapply(x, scale, 'scale', 1) ;

    idx_s_z = (ret.scale == 0) ;
    if (any(idx_s_z))
        warning ('Scale vaules == 0 detected. Changed to 1.')
        ret.scale(idx_s_z) = 1
    end

    n = size (x, 1) ;

    ret.x = (x - (ones (n, 1) * ret.center)) ./ (ones (n, 1) * ret.scale) ;
