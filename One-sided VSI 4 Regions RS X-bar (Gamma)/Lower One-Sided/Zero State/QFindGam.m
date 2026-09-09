% Create a fucntion to find the probability matrix Q %
% Lower One-sided VSI 4 Regions Run Sum X-bar Chart for Gamma %
% Zero State % 

function [Q, R1, u] = QFindGam(n, K, delta, S1, S2, S3, S4, a_shape, b_scale)

% Find transient-probability matrix Q 

B = zeros(4, 1) ; 
B(1, 1) = -S1 ; 
B(2, 1) = -S2 ; 
B(3, 1) = -S3 ; 
B(4, 1) = -S4 ; 

% Probability for each region 

mu = a_shape * b_scale ; % mean for xbar
var = a_shape * b_scale ^ 2 ; % variance for xbar

% Lower control limits 

LCL0 = max(0, mu - delta * sqrt(var)) ; 
LCL1 = max(0, mu - 1 * K * sqrt (var/n) - delta * sqrt(var)) ;
LCL2 = max(0, mu - 2 * K * sqrt (var/n) - delta * sqrt(var)) ;
LCL3 = max(0, mu - 3 * K * sqrt (var/n) - delta * sqrt(var)) ;

% Probability for each region 

p_3 = gamcdf(LCL2, n * a_shape, b_scale/n) - gamcdf(LCL3, n * a_shape, b_scale/n) ; 
p_2 = gamcdf(LCL1, n * a_shape, b_scale/n) - gamcdf(LCL2, n * a_shape, b_scale/n) ;
p_1 = gamcdf(LCL0, n * a_shape, b_scale/n) - gamcdf(LCL1, n * a_shape, b_scale/n) ;
p0 = 1 - gamcdf(LCL0, n * a_shape, b_scale/n) ; 

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
    for j = 5:columns
        if mod(j, 2) == 1
            R(a, j) = R(a, 1) + B(j - i, 1) ; 
            i = i + 1 ;
        else
            R(a, j) = 0 ; 
        end
    end
    for j = 3:columns
        if mod(j, 2) == 1
            if R(a, j) > -S4 
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
        if R1(t, 2) > R1(h, 2) 
            temp = R1(t, 2) ; 
            R1(t, 2) = R1(h, 2) ; 
            R1(h, 2) = temp ; 
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
Pr(1, 2) = p_1 ; 
Pr(1, 3) = p_2 ; 
Pr(1, 4) = p_3 ; 

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