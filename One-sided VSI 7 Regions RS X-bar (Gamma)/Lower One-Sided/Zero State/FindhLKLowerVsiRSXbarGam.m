% Create function to find hL and K based on ATS %
% Lower One-sided VSI 7 Regions Run Sum X-bar Chart for Gamma %
% Zero State %

function [hL, K] = FindhLKLowerVsiRSXbarGam(ATS0_target, ASI0_target, n, hS, G, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale)

    % Initial guess
    x0 = [0.8, 1.0];

    % Bounds for hL and K
    lb = [1.0, 0.6];  % Lower bounds for hL and K
    ub = [1.5, 2.5]; % Upper bounds for hL and K

    % Objective function with penalty for deviation from tau0
    obj_fun = @(x) constraints(x, ATS0_target, ASI0_target, n, hS, G, S1, S2, S3, S4, S5, S6,S7, a_shape, b_scale) ; 

    % Tolerance
    tol = 1e-15;

    % Options for lsqnonlin
    options = optimoptions('lsqnonlin', 'TolFun', tol, 'TolX', tol, 'MaxIter', 1000, 'Display', 'off');

    % Call to lsqnonlin to minimize the objective function with bounds
    [x, ~, ~] = lsqnonlin(obj_fun, x0, lb, ub, options);

    hL = x(1);
    K = x(2);

end