% deluppgiftB_basnivå.m
% Basnivå – deluppgift b) för Projekt B (vattenkran)
% Hitta K0 så att y'(L) = target_slope med egen RK4 och sekantmetod

% Parametrar och startvärden
K1           = 0.2;                     % Konstant K1
Y0           = [0.1; tan(deg2rad(46))]; % [y(0); y'(0)]
x0           = 0;                       % Startpunkt för x
L            = 0.5;                     % Slutpunkt för x
h            = 1e-3;                    % Steglängd för RK4
target_slope = -0.51;                   % Målvärde y'(L)

% Sekantmetod för att hitta K0
K0_prev = 5;  % Första gissning
K0_curr = 15; % Andra gissning

% Beräkna funktionsvärden f(K0) = y'(L;K0) - target_slope
f_prev = compute_end_slope(K0_prev, x0, Y0, h, L, K1) - target_slope;
f_curr = compute_end_slope(K0_curr, x0, Y0, h, L, K1) - target_slope;

tol_K   = 1e-5;  % Tolerans för K0
maxIter = 50;    % Max antal iterationer
iter    = 0;

while abs(K0_curr - K0_prev) > tol_K && iter < maxIter
    % Sekantuppdatering
    K0_new = K0_curr - f_curr*(K0_curr - K0_prev)/(f_curr - f_prev);
    % Flytta gissningar framåt
    K0_prev = K0_curr;
    f_prev  = f_curr;
    K0_curr = K0_new;
    f_curr  = compute_end_slope(K0_curr, x0, Y0, h, L, K1) - target_slope;
    iter = iter + 1;
end

K0_solution = K0_curr;  % Lösning för K0

% Metodfelsskattning med halverad steglängd
y_work = compute_end_slope(K0_solution, x0, Y0, h,    L, K1);
y_ref  = compute_end_slope(K0_solution, x0, Y0, h/2,  L, K1);
method_err = abs(y_ref - y_work);

% Visa resultat
fprintf('Löst K0          = %.6f\n', K0_solution);
fprintf('Antal iterationer= %d\n', iter);
fprintf('y''(L) med h     = %.6f\n', y_work);
fprintf('y''(L) med h/2   = %.6f\n', y_ref);
fprintf('Metodfel        = %.2e\n', method_err);

function s_end = compute_end_slope(K0, x0, Y0, h, L, K1)
    % Integrerar ODE med RK4 och returnerar y'(L)
    f = @(x,Y) ode_system(x, Y, K0, K1);
    [~, Y] = rk4_system(f, x0, Y0, h, L);
    s_end = Y(2,end);
end
