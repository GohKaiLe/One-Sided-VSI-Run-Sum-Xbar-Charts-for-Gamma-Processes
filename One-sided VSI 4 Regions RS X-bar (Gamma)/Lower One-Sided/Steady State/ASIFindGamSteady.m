% Create a fucntion to find ASI %
% Lower One-sided VSI 4 Regions Run Sum X-bar Chart for Gamma %
% Steady State % 

function [ASI] = ASIFindGamSteady(n, K, delta, hS, hL, G, S1, S2, S3, S4, a_shape, b_scale)

[Q, R1, u] = QFindGamSteady(n, K, delta, S1, S2, S3, S4, a_shape, b_scale) ; 

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

one = ones(u, 1) ; 

s0 = ((pinv(I - Q')) * s ) / (one' * (pinv(I - Q')) * s) ;

alpha = zeros(u, 1) ; 

for i = 1:u 
    alpha(i) = (s0(i) * h(i)) / (s0' * h) ; 
end

ASI = (s0' * (pinv(I - Q)) * h) / (s0' * (pinv(I - Q)) * ones(u, 1)) ; 


end