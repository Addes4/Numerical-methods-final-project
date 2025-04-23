function [slope_end, height_max] = compute_metrics(K0, s0, x0, Y0, h, L, K1)
% COMPUTE_METRICS Returnerar y'(L) och max y över [x0,L].
    Y0_loc = [Y0(1); s0];
    f = @(x,Y) crane_ode(x, Y, K0, K1);
    [~, Y] = rk4_system(f, x0, Y0_loc, h, L);
    slope_end  = Y(2,end);
    height_max = max(Y(1,:));
end
