function [ ret ] = p_genapply ( x, arg, arg_name, def )

    p = size (x, 2) ;
   
    if (isempty (arg)) 
        ret = ones (1, p) * def ;
    elseif (isnumeric(arg))
        if (length(arg) == 1) 
            ret = ones (p, 1) * arg ;
        elseif (length (arg) == p)
            ret = arg ;
        else
            error ('Argument arg has wrong length. Expected either 1 or p.') ;
        end
    elseif (ischar (arg))
        f = str2func (arg) ;

        ret = f (x) ;

        if (length (ret) ~= p)
            error ('Function arg returned vector of wrong length.')
        end
    else
        error ('Argument arg is of unknown type') ;
    end
