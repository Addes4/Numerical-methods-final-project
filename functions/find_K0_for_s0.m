function K0_out = find_K0_for_s0(s0, target_slope, x0, Y0, h, L, K1)
    % initiera gissningar
    K0_prev = 5;
    K0_curr = 15;
    tol_K0  = 1e-4;
    maxIter = 50;
    [s_prev, ~] = compute_metrics(K0_prev, s0, x0, Y0, h, L, K1);
    f_prev = s_prev - target_slope;
    [s_curr, ~] = compute_metrics(K0_curr, s0, x0, Y0, h, L, K1);
    f_curr = s_curr - target_slope;
    iter = 0;
    % sekantloop
    while abs(K0_curr - K0_prev) > tol_K0 && iter < maxIter
        K0_new   = K0_curr - f_curr*(K0_curr - K0_prev)/(f_curr - f_prev);
        K0_prev  = K0_curr;  f_prev = f_curr;
        K0_curr  = K0_new;
        [s_curr, ~] = compute_metrics(K0_curr, s0, x0, Y0, h, L, K1);
        f_curr  = s_curr - target_slope;
        iter = iter + 1;
    end
    K0_out = K0_curr;
end
