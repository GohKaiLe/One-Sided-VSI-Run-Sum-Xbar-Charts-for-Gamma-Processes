% Create a function to find SDTS %
% Lower One-sided VSI 4 Regions Run Sum X-bar Chart for Gamma %
% Zero State % 

function [SDTS] = SDTSFindGam(n, K, delta, hS, hL, G, S1, S2, S3, S4, a_shape, b_scale)

[Q, R1, u] = QFindGam(n, K, delta, S1, S2, S3, S4, a_shape, b_scale) ; 

s = zeros(u, 1) ; 
s(1, 1) = 1 ; 
I = eye(u) ; 

H = S4 ; 
h = ones(u, 1) ;  
temp3 = R1(:, 1) ; % extract all possible lower cumulativen scores

for i = 1:u
    if (-H < temp3(i)) && (temp3(i) <= -H/G)
        h(i) = hS ;
    else
        h(i) = hL ;
    end
end

E = diag(h) ; 

v1 = s' * (inv(I - Q)) * E * (2 * (inv(I - Q)) - I) * h ; 
v2 = (s' * (inv(I - Q)) * h)^2 ; 

SDTS = sqrt(v1 - v2) ; 


end