
%parametrar
K1= 0.2; % Konstant K1
Y0 = [0.1; tan(deg2rad(46))]; % [y(0); y'(0)]
x0 = 0; 
L = 0.5;
h = 1e-3;

% Mål: justera K0 så att y'(L) = target_slope
target_slope = -0.51;
% För root-finding via sekantmetod, två startgissningar på K0
K0_prev = 5; % gissning 1
K0_curr = 15; % gissning 2
% Beräkna funktionsvärden med grov steglängd h
f_prev = compute_end_slope(K0_prev, x0, Y0, h, L, K1) - target_slope;
f_curr = compute_end_slope(K0_curr, x0, Y0, h, L, K1) - target_slope;

tol_K = 1e-5; %tolerans för K0
maxIter = 50;
iter = 0;

while abs(K0_curr - K0_prev) > tol_K && iter < maxIter
    % Sekantuppdatering
    K0_new   = K0_curr - f_curr*(K0_curr - K0_prev)/(f_curr - f_prev);
    % Uppdatera för nästa iteration
    K0_prev = K0_curr;    f_prev = f_curr;
    K0_curr = K0_new;
    f_curr = compute_end_slope(K0_curr, x0, Y0, h, L, K1) - target_slope;
    iter = iter + 1;
end

K0_solution = K0_curr;  % lösta K0 värdet

% Beräkna metodfel med halverad steglängd h/2
y_work = compute_end_slope(K0_solution, x0, Y0, h,    L, K1);
y_ref = compute_end_slope(K0_solution, x0, Y0, h/2,  L, K1);
method_err = abs(y_ref - y_work);

function s_end = compute_end_slope(K0, x0, Y0, h, L, K1)
    % Räknar ut y'(L) med RK4 för givna parametrar
    f = @(x,Y) crane_ode(x, Y, K0, K1);
    [~, Y] = rk4_system(f, x0, Y0, h, L);
    s_end  = Y(2,end);
end

function dYdx = crane_ode(x, Y, K0, K1)
    % ODE-system: Y = [y; y']
    y = Y(1);
    v = Y(2);
    dYdx = [
        v;
        -(K0 - K1*x)*y*(1 + v^2)^(3/2)
        ];
end

function [x, Y] = rk4_system(f, x0, Y0, h, L)
    % Fyra-stegs Runge–Kutta
    x = x0:h:L;
    N = numel(x);
    Y = zeros(2, N);
    Y(:,1)= Y0;
    for i = 1:N-1
        k1 = f(x(i), Y(:,i));
        k2 = f(x(i)+h/2, Y(:,i) + h*k1/2);
        k3 = f(x(i)+h/2, Y(:,i) + h*k2/2);
        k4 = f(x(i)+h, Y(:,i) + h*k3);
        Y(:,i+1)= Y(:,i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end

disp(['nya K0 = ' num2str(K0_solution,'%.6f')])
disp(['metodfelet = ' num2str(method_err,'%.2e')])
disp(['y''(0.5) (h)  = ' num2str(y_work,'%.6f')])
disp(['y''(0.5) (h/2)= ' num2str(y_ref,'%.6f')])
