% Create a fucntion to find ASI %
% Upper One-sided VSI 4-region RS X-bar Chart for the Gamma Distribution %

function [ASI] = ASIFindGam(n, K, delta, hS, hL, G, S1, S2, S3, S4, a_shape, b_scale)

[Q, R1, u] = QFindGam(n, K, delta, S1, S2, S3, S4, a_shape, b_scale) ; 

s = zeros(u, 1) ; 
s(1, 1) = 1 ; 
I = eye(u) ; 

H = S4 ; 
h = ones(u, 1) ;  
temp2 = R1(:, 1) ; % extract all possible upper cumulativen scores

for i = 1:u
    if (H/G <= temp2(i)) && (temp2(i) < H)
        h(i) = hS ;
    else
        h(i) = hL ;
    end 
end

ASI = (s' * (inv(I - Q)) * h) / (s' * (inv(I - Q)) * ones(u, 1)) ; 

end