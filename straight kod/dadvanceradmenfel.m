% Skript: deluppgift_d_basnivå.m
% Basnivå deluppgift d) – Använder MATLAB:s inbyggda funktioner för sekantmetoden (fzero)
% Använder extern fil rk4_system.m för RK4

% Parametrar och mål
x0             = 0;                      % Startpunkt x
L              = 0.5;                    % Slutpunkt x
Y0             = [0.1; tan(deg2rad(46))]; % [y(0); y'(0)]
K1             = 0.2;                    % Konstant K1
h              = 1e-5;                   % Steglängd för RK4

target_slope   = -0.51;                  % Målvärde y'(L)
target_height  = 0.255;                  % Målvärde max y

tol_K0         = 1e-5;                   % Tolerans för K0
tol_s0         = 1e-5;                   % Tolerans för s0

% 1) Justera startlutning s0 med fzero för att uppnå maxhöjd
s0_initial_guess1 = tan(deg2rad(40));
s0_initial_guess2 = tan(deg2rad(50));

% Definiera funktion för y_max - target_height
fun_outer = @(s0) compute_height_error(s0, h, x0, L, Y0, K1, target_slope, target_height, tol_K0);
options_outer = optimset('TolX', tol_s0);
s0_solution = fzero(fun_outer, [s0_initial_guess1, s0_initial_guess2], options_outer);

% 2) Bestäm K0 för det framräknade s0
K0_solution = find_K0_for_s0(s0_solution, h, x0, L, Y0, K1, target_slope, tol_K0);

% 3) Metodfelsskattning med h och h/2
[slope_h, height_h]   = compute_metrics(K0_solution, s0_solution, h, x0, L, Y0, K1);
[slope_h2, height_h2] = compute_metrics(K0_solution, s0_solution, h/2, x0, L, Y0, K1);
method_err_slope      = abs(slope_h2 - slope_h);
method_err_height     = abs(height_h2 - height_h);

% Visning av resultat
disp('Basnivå deluppgift d)')
disp(['  Startlutning s0   = ' num2str(s0_solution, '%.6f')])
disp(['  Parameter K0      = ' num2str(K0_solution, '%.6f')])
disp(['  Slutlutning y''   = ' num2str(slope_h, '%.6f')])
disp(['  Maxhöjd y_max     = ' num2str(height_h, '%.6f') ' m'])
disp(['  Metodfel lutning = ' num2str(method_err_slope, '%.2e')])
disp(['  Metodfel höjd    = ' num2str(method_err_height, '%.2e')])
disp(['  Startlutning s0 i grader = ' num2str(rad2deg(atan(s0_solution)), '%.6f')])

% Lokala hjälpfunktioner
function [s_end, y_max] = compute_metrics(K0, s0, h_step, x0, L, Y0, K1)
    odeFun = @(x, Y) [Y(2); -(K0 - K1*x)*Y(1)*(1+Y(2)^2)^(3/2)];
    [~, Y] = rk4_system(odeFun, x0, [Y0(1); s0], h_step, L);
    s_end = Y(2, end);
    y_max = max(Y(1, :));
end

function K0_out = find_K0_for_s0(s0, h_step, x0, L, Y0, K1, tgt_slope, tol_K0)
    fun = @(K0) compute_metrics(K0, s0, h_step, x0, L, Y0, K1) - tgt_slope;
    options = optimset('TolX', tol_K0);
    K0_out = fzero(fun, [5, 15], options);
end

function error = compute_height_error(s0, h_step, x0, L, Y0, K1, tgt_slope, tgt_height, tol_K0)
    K0 = find_K0_for_s0(s0, h_step, x0, L, Y0, K1, tgt_slope, tol_K0);
    [~, y_max] = compute_metrics(K0, s0, h_step, x0, L, Y0, K1);
    error = y_max - tgt_height;
end
