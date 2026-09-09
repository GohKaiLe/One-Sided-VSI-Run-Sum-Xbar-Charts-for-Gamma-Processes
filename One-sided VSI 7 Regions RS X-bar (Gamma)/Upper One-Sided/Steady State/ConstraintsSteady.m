% Constraint for AATS optimisation %
% Upper One-sided VSI 7 Regions Run Sum X-bar Chart for Gamma %
% Steady State %

function ceq = ConstraintsSteady(x, AATS0_target, ASI0_target, n, hS, G, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale)

    hL = x(1);
    K = x(2);

    AATS0 = AATSFindGam(n, K, 0, hS, hL, G, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale) ;
    ASI0 = ASIFindGamSteady(n, K, 0, hS, hL, G, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale) ;

    % Equality constraints
    ceq(1) = AATS0 - AATS0_target;
    ceq(2) = ASI0 - ASI0_target;

end