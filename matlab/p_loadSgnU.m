function [ x ] = p_loadSgnU ( x )
    [val.max, idx.max] = max (abs (x)) ;
    %idx.max <- apply(abs(x), 2, which.max) ;
    sgn = sign (p_MatSSperCol (x, idx.max)) ;
    %sgn = sign(x[cbind(idx.max, 1:ncol(x))]) ;
    x = x * diag(sgn) ;

end



