% Constraint for ATS optimisation %
% Lower One-sided VSI 7 Regions Run Sum X-bar Chart for Gamma %
% Zero State %

function ceq = constraints(x, ATS0_target, ASI0_target, n, hS, G, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale)

    hL = x(1);
    K = x(2);

    ATS0 = ATSFindGam(n, K, 0, hS, hL, G, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale) ;
    ASI0 = ASIFindGam(n, K, 0, hS, hL, G, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale) ;

    % Check for Inf/NaN
    if isnan(ATS0) || isinf(ATS0) || isnan(ASI0) || isinf(ASI0)
        ceq = [1e10; 1e10];  % Penalize invalid values
    else
        ceq = [ATS0 - ATS0_target; ASI0 - ASI0_target];
    end
end