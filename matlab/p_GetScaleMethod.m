function [ nScale ] = p_GetScaleMethod (szScale)
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here

nScale = -1 ;
if (ischar (szScale))
    if (strcmp (szScale, 'sd'))
       nScale = 0 ;
    elseif (strcmp (szScale, 'mad'))
       nScale = 1 ;
    elseif (strcmp (szScale, 'qn'))
       nScale = 2 ;
    end
elseif (isnumeric (szScale))
    if (isfinite (szScale))
        if (szScale >= 0 && szScale <= 2)
            nScale = szScale ;
        end
    end
end        

if (nScale < 0)
    error ('Unknown  scale provided..') ;
end
