% Create a fucntion to find the probability matrix Q %
% Upper One-sided VSI 7 Regions RS X-bar Chart for Gamma %
% Zero State %

function [Q, R1, u] = QFindGam(n, K, delta, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale)

a1 = S1 ; % Score region 1 (above mu)
a2 = S2 ; % Score region 2 (above mu)
a3 = S3 ; % Score region 3 (above mu)
a4 = S4 ; % Score region 4 (above mu)
a5 = S5 ; % Score region 5 (above mu)
a6 = S6 ; % Score region 6 (above mu)
a7 = S7 ; % Score region 7 (above mu)

A = zeros(7, 1) ;
A(1, 1) = a1 ;
A(2, 1) = a2 ;
A(3, 1) = a3 ;
A(4, 1) = a4 ;
A(5, 1) = a5 ; 
A(6, 1) = a6 ; 
A(7, 1) = a7 ; 

% Probability for each region 

mu = a_shape * b_scale ; % mean for xbar
var = a_shape * b_scale ^ 2 ; % variance for xbar

% Upper control limits 

UCL0 = mu - delta * sqrt(var) ; 
UCL1 = mu + 0.5 * K * sqrt(var/n) - delta * sqrt(var) ; 
UCL2 = mu + 1.0 * K * sqrt(var/n) - delta * sqrt(var) ; 
UCL3 = mu + 1.5 * K * sqrt(var/n) - delta * sqrt(var) ; 
UCL4 = mu + 2.0 * K * sqrt(var/n) - delta * sqrt(var) ; 
UCL5 = mu + 2.5 * K * sqrt(var/n) - delta * sqrt(var) ; 
UCL6 = mu + 3.0 * K * sqrt(var/n) - delta * sqrt(var) ; 

% Probability for each region

p6 = gamcdf(UCL6, n * a_shape, b_scale/n) - gamcdf(UCL5, n * a_shape, b_scale/n) ; 
p5 = gamcdf(UCL5, n * a_shape, b_scale/n) - gamcdf(UCL4, n * a_shape, b_scale/n) ; 
p4 = gamcdf(UCL4, n * a_shape, b_scale/n) - gamcdf(UCL3, n * a_shape, b_scale/n) ; 
p3 = gamcdf(UCL3, n * a_shape, b_scale/n) - gamcdf(UCL2, n * a_shape, b_scale/n) ; 
p2 = gamcdf(UCL2, n * a_shape, b_scale/n) - gamcdf(UCL1, n * a_shape, b_scale/n) ; 
p1 = gamcdf(UCL1, n * a_shape, b_scale/n) - gamcdf(UCL0, n * a_shape, b_scale/n) ; 
p0 = gamcdf(UCL0, n * a_shape, b_scale/n) ;

rows = 50 ; 
columns = 18 ; 
R = zeros(rows, columns) ; 
S = zeros(rows, 2) ; 
S(1, 1) = 0 ; 
S(1, 2) = 0 ; 

u = 1 ; 

for a = 1:rows
    R(a, 1) = S(a, 1) ; 
    R(a, 2) = S(a, 2) ; 
    for j = 3:4 
        R(a, j) = 0 ;
    end
    i = 4 ;
    for j = 5:columns
        if mod(j, 2) == 1 
            R(a, j) = R(a, 1) + A(j - i, 1) ; 
            i = i + 1 ;
        else
            R(a, j) = 0 ;
        end
    end
    for j = 3:columns 
        if mod(j, 2) == 1
            if R(a, j) < S7
                ww = 0 ;
                for w = 1:rows
                    if R(a, j) ~= S(w, 1) 
                        ww = ww + 1 ; 
                    end
                end
                if ww == rows 
                    S(u + 1, 1) = R(a, j) ;
                    S(u + 1, 2) = 0 ;
                    u = u + 1 ;
                    ww = 0 ;
                end
            end
        end
    end
end

M = zeros(u, rows) ; 

for i = 1:u 
    M(i, i) = 1 ; 
end

R1 = M * R ; 
temp1 = zeros(1, columns) ; 

for t = 2:(u - 1)
    for h = (t + 1):u
        if R1(t, 1) > R1(h, 1) 
            temp = R1(t, 1) ; 
            R1(t, 1) = R1(h, 1) ; 
            R1(h, 1) = temp ; 
            for i = 1:columns
                if i ~= 1
                    temp1(1, i) = R1(t, i) ; 
                    R1(t, i) = R1(h, i) ; 
                    R1(h, i) = temp1(1, i) ; 
                end
            end
        end
    end
end

R2 = zeros(u, 9) ; 

for b = 1:u 
    c = 1 ;
    for a = 3:2:columns 
        d = 0 ;
        for row = 1:u
            if R1(b, a) == R1(row, 1) && R1(b, a + 1) == R1(row, 2) 
                R2(b, a - c) = row ; 
                d = 1 ;
            end
        end
        if d ~= 1 
            R2(b, a - c) = u + 1 ;
        end
        c = c + 1 ; 
    end
end

for a = 1:u
    R2(a, 1) = a ; 
end

Q = zeros(u, u) ; 
Pr = zeros(1, 7) ; 
Pr(1, 1) = p0 ; 
Pr(1, 2) = p1 ; 
Pr(1, 3) = p2 ; 
Pr(1, 4) = p3 ; 
Pr(1, 5) = p4 ; 
Pr(1, 6) = p5 ; 
Pr(1, 7) = p6 ; 

for a = 1:u
    for p = 1:u
        prob = 0 ;
        for b = 2:8
            if R2(a, b) == p
                prob = prob + Pr(1, b - 1) ; 
            end
        end
        Q(a, p) = prob ; 
    end
end

end