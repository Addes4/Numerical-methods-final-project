% Skript: deluppgift_d_basnivå.m
% Basnivå deluppgift d) – Egen RK4 och sekantmetod med metodfelsskattning
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
maxIter_K0     = 50;                     % Max iterationer för K0-sökning
tol_s0         = 1e-5;                   % Tolerans för s0
maxIter_s0     = 50;                     % Max iterationer för s0-sökning

% 1) Justera startlutning s0 med sekantmetod för maxhöjd
s0_prev = tan(deg2rad(40));
s0_curr = tan(deg2rad(50));
it = 0;
while abs(s0_curr - s0_prev) > tol_s0 && it < maxIter_s0
    % Beräkna maxhöjd för två gissningar på s0
    K0_prev = find_K0_for_s0(s0_prev, h, x0, L, Y0, K1, target_slope, tol_K0, maxIter_K0);
    [~, h_prev] = compute_metrics(K0_prev, s0_prev, h, x0, L, Y0, K1);
    K0_curr = find_K0_for_s0(s0_curr, h, x0, L, Y0, K1, target_slope, tol_K0, maxIter_K0);
    [~, h_curr] = compute_metrics(K0_curr, s0_curr, h, x0, L, Y0, K1);
    % Sekantsteg för s0
    f_prev = h_prev - target_height;
    f_curr = h_curr - target_height;
    s_new  = s0_curr - f_curr*(s0_curr - s0_prev)/(f_curr - f_prev);
    s0_prev = s0_curr;
    s0_curr = s_new;
    it = it + 1;
end
s0_solution = s0_curr;

% 2) Bestäm K0 för det framräknade s0
K0_solution = find_K0_for_s0(s0_solution, h, x0, L, Y0, K1, target_slope, tol_K0, maxIter_K0);

% 3) Metodfelsskattning med h och h/2
[slope_h, height_h]   = compute_metrics(K0_solution, s0_solution,    h,  x0, L, Y0, K1);
[slope_h2,height_h2]   = compute_metrics(K0_solution, s0_solution,   h/2, x0, L, Y0, K1);
method_err_slope      = abs(slope_h2   - slope_h);
method_err_height     = abs(height_h2  - height_h);

% Lokala hjälpfunktioner (fortfarande i samma fil)

function [s_end, y_max] = compute_metrics(K0, s0, h_step, x0, L, Y0, K1)
    % Beräknar slutlutning och maxhöjd med RK4 via extern rk4_system
    odeFun = @(x, Y) [Y(2); -(K0 - K1*x)*Y(1)*(1+Y(2)^2)^(3/2)];
    [~, Y] = rk4_system(odeFun, x0, [Y0(1); s0], h_step, L);
    s_end = Y(2,end);
    y_max  = max(Y(1,:));
end

function K0_out = find_K0_for_s0(s0, h_step, x0, L, Y0, K1, tgt_slope, tol_K0, maxIter_K0)
    % Sekantmetod för K0 så att y'(L)=tgt_slope
    K0_prev = 5;    K0_curr = 15;    iter=0;
    f_prev  = compute_metrics(K0_prev, s0, h_step, x0, L, Y0, K1) - tgt_slope;
    f_curr  = compute_metrics(K0_curr, s0, h_step, x0, L, Y0, K1) - tgt_slope;
    while abs(K0_curr-K0_prev)>tol_K0 && iter<maxIter_K0
        K0_new = K0_curr - f_curr*(K0_curr-K0_prev)/(f_curr-f_prev);
        K0_prev=K0_curr; f_prev=f_curr;
        K0_curr=K0_new;
        f_curr = compute_metrics(K0_curr, s0, h_step, x0, L, Y0, K1) - tgt_slope;
        iter=iter+1;
    end
    K0_out = K0_curr;
end

disp('Basnivå deluppgift d)')
disp(['  Startlutning s0   = ' num2str(s0_solution,   '%.6f')])
disp(['  Parameter K0      = ' num2str(K0_solution,   '%.6f')])
disp(['  Slutlutning y''   = ' num2str(slope_h,        '%.6f')])
disp(['  Maxhöjd y_max     = ' num2str(height_h,       '%.6f') ' m'])
disp(['  Metodfel lutning = ' num2str(method_err_slope,'%.2e')])
disp(['  Metodfel höjd    = ' num2str(method_err_height,'%.2e')])
disp(['  Startlutning s0   = ' num2str(rad2deg(s0_solution),   '%.6f')])
