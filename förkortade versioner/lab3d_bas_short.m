% Skript: deluppgift_d_basnivå.m
% Basnivå deluppgift d) – Egen RK4 och dubbel sekantmetod med metodfelsskattning
% Använder externa filer:
%   • rk4_system.m
%   • ode_system.m
%   • compute_metrics.m
%   • find_K0_for_s0.m


% Parametrar och mål
x0            = 0;                          % Startpunkt x
L             = 0.5;                        % Slutpunkt x
Y0            = [0.1; tan(deg2rad(46))];    % [y(0); y'(0)]
K1            = 0.2;                        % Konstant K1
h             = 1e-5;                       % Steglängd för RK4

target_slope  = -0.51;                      % Målvärde y'(L)
target_height = 0.255;                      % Målvärde max y

tol_s0        = 1e-5;                       % Tolerans för s0
maxIter_s0    = 50;                         % Max iterationer för s0
tol_K0        = 1e-5;                       % Tolerans för K0
maxIter_K0    = 50;                         % Max iterationer för K0

% 1) Yttre sekant: justera s0 för att få rätt maxhöjd
s0_prev = tan(deg2rad(40));
s0_curr = tan(deg2rad(50));
it_s     = 0;

while abs(s0_curr - s0_prev) > tol_s0 && it_s < maxIter_s0
    % Beräkna maxhöjd med nuvarande gissningar
    K0_a   = find_K0_for_s0(s0_prev,  target_slope, x0, Y0, h, L, K1);
    [~, h_a] = compute_metrics( K0_a,   s0_prev,   x0, Y0, h, L, K1);
    K0_b   = find_K0_for_s0(s0_curr,  target_slope, x0, Y0, h, L, K1);
    [~, h_b] = compute_metrics( K0_b,   s0_curr,   x0, Y0, h, L, K1);

    % Sekantsteg för s0
    f_a = h_a - target_height;
    f_b = h_b - target_height;
    s_new = s0_curr - f_b*(s0_curr - s0_prev)/(f_b - f_a);

    % Uppdatera
    s0_prev = s0_curr;
    s0_curr = s_new;
    it_s    = it_s + 1;
end

s0_solution = s0_curr;

% 2) Inre sekant: bestäm K0 för funnet s0
K0_solution = find_K0_for_s0(...
    s0_solution, target_slope, x0, Y0, h, L, K1);

% 3) Metodfelsskattning med h och h/2
[slope_h,  height_h ] = compute_metrics( ...
    K0_solution, s0_solution, x0, Y0, h,  L, K1);
[slope_h2, height_h2] = compute_metrics( ...
    K0_solution, s0_solution, x0, Y0, h/2, L, K1);

err_slope  = abs(slope_h2  - slope_h);
err_height = abs(height_h2 - height_h);

% 4) Visa resultat
disp('Basnivå – deluppgift d)');
fprintf('  Startlutning s0   = %.6f rad (%.4f°)\n', ...
    s0_solution, rad2deg(s0_solution));
fprintf('  Parameter K0      = %.6f\n', K0_solution);
fprintf('  Slutlutning y''   = %.6f\n', slope_h);
fprintf('  Maxhöjd y_max     = %.6f m\n', height_h);
fprintf('  Metodfel lutning = %.2e\n', err_slope);
fprintf('  Metodfel höjd    = %.2e\n', err_height);
