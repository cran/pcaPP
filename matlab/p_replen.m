function [ y ] = p_replen ( x, n )

    six = size (x) ;
    if (length (six) > 2 || ~any (six == 1))
        error ('x must be a vector') ;
    end
    if (six(1) == 1)
        x = x' ;
    end    
       
    lx = length (x) ;
    
    if (lx >= n)
        y = x (1:n) ;
    else
        y = zeros (1, n) ;
        for i = 1:floor (n / lx)
            y ( ((i-1)*lx + 1) : (i*lx)) = x ;
        end
        
        di = n - i * lx  ;
        
        y (i*lx +1 : n) = x (1:di) ;

    end