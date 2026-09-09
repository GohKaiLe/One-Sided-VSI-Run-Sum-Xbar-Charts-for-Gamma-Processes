% Create a function to find K and hL %
% Upper One-sided VSI 4 Regions Run Sum X-bar Chart for Gamma %

function [hL, K] = FindhLKUpperVsiRSXbarGam(ATS0_target, ASI0_target, n, hS, G, S1, S2, S3, S4, a_shape, b_scale)

    % Initial guess
    x0 = [1.2, 2.0];

    % Bounds for hL and K
    lb = [1.0, 0.6];  % Lower bounds for hL and K
    ub = [1.5, 3.0]; % Upper bounds for hL and K

    % Objective function with penalty for deviation from tau0
    obj_fun = @(x) constraints(x, ATS0_target, ASI0_target, n, hS, G, S1, S2, S3, S4, a_shape, b_scale) ; 

    % Tolerance
    tol = 1e-15;

    % Options for lsqnonlin
    options = optimoptions('lsqnonlin', 'TolFun', tol, 'TolX', tol, 'MaxIter', 2000, 'Display', 'off');

    % Call to lsqnonlin to minimize the objective function with bounds
    [x, ~, ~] = lsqnonlin(obj_fun, x0, lb, ub, options);

    hL = x(1);
    K = x(2);

end