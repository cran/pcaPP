function [ y ] = p_MatSSperCol ( x, id )

    p = size (x, 2) ;
    y = zeros (1, p) ;

    if( length (id) ~= p)
       error ('id has to be of length size (x, 2).') 
    end

    for (i = 1:p)
        y(i) = x (id(i), i) ;
    end

end

