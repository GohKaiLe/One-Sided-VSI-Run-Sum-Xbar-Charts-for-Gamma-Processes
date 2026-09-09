% Create a fucntion to find ASI %
% Upper One-sided VSI 7 Regions Run Sum X-bar Chart for Gamma %
% Steady State %

function [ASI] = ASIFindGamSteady(n, K, delta, hS, hL, G, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale)

[Q, R1, u] = QFindGamSteady(n, K, delta, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale) ; 

s = zeros(u, 1) ; 
s(1, 1) = 1 ; 
I = eye(u) ; 

H = S7 ; 
h = ones(u, 1) ;  
temp2 = R1(:, 1) ; % extract all possible upper cumulativen scores

for i = 1:u
    if (H/G <= temp2(i)) && (temp2(i) < H)
        h(i) = hS ;
    else
        h(i) = hL ;
    end 
end

one = ones(u, 1) ; 

s0 = ((inv(I - Q')) * s ) / (one' * (inv(I - Q')) * s) ;

ASI = (s0' * (inv(I - Q)) * h) / (s0' * (inv(I - Q)) * ones(u, 1)) ; 

end