% Create a function to find EAATS %
% Upper One-sided VSI 7 Regions Run Sum X-bar Chart for Gamma %
% Steady State %

function [EAATS] = EAATSFindUpperGam(n, K, deltamin, deltamax, hS, hL, G, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale)

h = 9 ; 

[x, w] = lgwt_table(h) ; 

EAATSsum = 0 ; 

for i = 1:h 
    xi = x(i) ; 
    wi = w(i) ; 
    d = ((deltamax - deltamin)/2) * xi + ((deltamax + deltamin)/2) ; 
    AATS = AATSFindGam(n, K, d, hS, hL, G, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale) ; 
    EAATS = (1/(deltamax - deltamin)) * AATS * wi ; 
    EAATSsum = EAATSsum + EAATS ; 
end

EAATS = ((deltamax - deltamin)/2) * EAATSsum ; 

end