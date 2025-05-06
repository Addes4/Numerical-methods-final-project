
% parameter och mål
x0 = 0;
L = 0.5;
Y0 = [0.1; tan(deg2rad(46))]; % [y(0); y'(0)]
K1  = 0.2;
h = 1e-3;

target_slope  = -0.51;
target_height = 0.255;

tol_s0 = 1e-3;
maxIter_s0 = 50;
tol_K0 = 1e-3;
maxIter_K0 = 50;

% yttre sekant
s0_prev = tan(deg2rad(40));
s0_curr = tan(deg2rad(50));
it_s = 0;

while abs(s0_curr - s0_prev) > tol_s0 && it_s < maxIter_s0
    %  maxhöjd med nuvarande gissningar
    K0_a   = find_K0_for_s0(s0_prev,  target_slope, x0, Y0, h, L, K1);
    [~, h_a] = compute_metrics( K0_a,   s0_prev,   x0, Y0, h, L, K1);
    K0_b   = find_K0_for_s0(s0_curr,  target_slope, x0, Y0, h, L, K1);
    [~, h_b] = compute_metrics( K0_b,   s0_curr,   x0, Y0, h, L, K1);

    % sekantsteg för s0
    f_a = h_a - target_height;
    f_b = h_b - target_height;
    s_new = s0_curr - f_b*(s0_curr - s0_prev)/(f_b - f_a);

    % upd
    s0_prev = s0_curr;
    s0_curr = s_new;
    it_s    = it_s + 1;
end

s0_solution = s0_curr;

% inre sekant
K0_solution = find_K0_for_s0(s0_solution, target_slope, x0, Y0, h, L, K1);

% metodfelsskattning
[slope_h,  height_h ] = compute_metrics(K0_solution, s0_solution, x0, Y0, h,  L, K1);
[slope_h2, height_h2] = compute_metrics(K0_solution, s0_solution, x0, Y0, h/2, L, K1);

err_slope  = abs(slope_h2  - slope_h);
err_height = abs(height_h2 - height_h);

fprintf('startlutning s0 = %.6f rad (%.4f grader)\n', s0_solution, rad2deg(s0_solution));
fprintf('parameter K0 = %.6f\n', K0_solution);
fprintf('slutlutning y'' = %.6f\n', slope_h);
fprintf('maxhöjd y_max = %.6f m\n', height_h);
fprintf('metodfel lutning = %.2e\n', err_slope);
fprintf('metodfel höjd = %.2e\n', err_height);
