% Create a fucntion to find the probability matrix Q %
% Upper One-sided VSI 4-region RS X-bar Chart for the Gamma Distribution % 

function [Q, R1, u] = QFindGam(n, K, delta, S1, S2, S3, S4, a_shape, b_scale)

% Find transient-probability matrix Q 

A = zeros(4, 1) ;
A(1, 1) = S1 ;
A(2, 1) = S2 ;
A(3, 1) = S3 ;
A(4, 1) = S4 ;

% Probability for each region 

mu = a_shape * b_scale ; % mean for xbar
var = a_shape * b_scale ^ 2 ; % variance for xbar

UCL0 = mu - delta * sqrt(var) ; 
UCL1 = mu + K * sqrt (var/n) - delta * sqrt(var) ;
UCL2 = mu + 2 * K * sqrt (var/n) - delta * sqrt(var) ;
UCL3 = mu + 3 * K * sqrt (var/n) - delta * sqrt(var) ;

p3 = gamcdf(UCL3, n * a_shape, b_scale / n) - gamcdf(UCL2, n * a_shape, b_scale / n) ;
p2 = gamcdf(UCL2, n * a_shape, b_scale / n) - gamcdf(UCL1, n * a_shape, b_scale / n) ;
p1 = gamcdf(UCL1, n * a_shape, b_scale / n) - gamcdf(UCL0, n * a_shape, b_scale / n) ;
p0 = gamcdf(UCL0, n * a_shape, b_scale / n) ; 

rows = 50 ; 
columns = 12 ; 
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
    for j =5:columns
        if mod(j, 2) == 1 
            R(a, j) = R(a, 1) + A(j - i, 1) ; 
            i = i + 1 ;
        else
            R(a, j) = 0 ;
        end
    end
    for j = 3:columns 
        if mod(j, 2) == 1
            if R(a, j) < S4
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

R2 = zeros(u, 6) ; 

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
Pr = zeros(1, 4) ; 
Pr(1, 1) = p0 ; 
Pr(1, 2) = p1 ; 
Pr(1, 3) = p2 ; 
Pr(1, 4) = p3 ; 

for a = 1:u
    for p = 1:u
        prob = 0 ;
        for b = 2:5
            if R2(a, b) == p
                prob = prob + Pr(1, b - 1) ; 
            end
        end
        Q(a, p) = prob ; 
    end
end

end