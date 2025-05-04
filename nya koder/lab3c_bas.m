%utan interpolation

%parametrar
K1 = 0.2; 
Y0= [0.1; tan(deg2rad(46))];   % [y(0); y'(0)]
x0= 0; 
L  = 0.5;
h = 1e-3; %steglängd för RK4

%bestäm K0 fr deluppgift b) m sekant
target_slope  = -0.51; %målvärde för y'(L)
K0_prev = 5; % gissning 1
K0_curr = 15; % gissning 2
f_prev = compute_slope(K0_prev, x0, Y0, h, L, K1) - target_slope;
f_curr = compute_slope(K0_curr, x0, Y0, h, L, K1) - target_slope;

tol_K = 1e-4; % tolerans för K0
maxIter = 50; % max iterationer
iter = 0;

while abs(K0_curr - K0_prev) > tol_K && iter < maxIter
    K0_new = K0_curr - f_curr * (K0_curr - K0_prev) / (f_curr - f_prev);
    K0_prev = K0_curr;  f_prev = f_curr;
    K0_curr = K0_new;
    f_curr = compute_slope(K0_curr, x0, Y0, h, L, K1) - target_slope;
    iter = iter + 1;
end
K0_solution = K0_curr; % lös för K0

%beräkna kranprofil och maxhöjd
% m steglängd h
f = @(x,Y) crane_ode(x, Y, K0_solution, K1);
[~, Yh ] = rk4_system(f, x0, Y0, h,    L);
max_yh  = max(Yh(1,:)); %maxhöjd för h

%m steglängd h/2 för felet
h2 = h/2;
[~, Yh2] = rk4_system(f, x0, Y0, h2,   L);
max_yh2 = max(Yh2(1,:)); % maxhöjd för h/2

method_err = abs(max_yh2 - max_yh);

function slope_end = compute_slope(K0, x0, Y0, h, L, K1)
    %räknar y'(L) med RK4 för given K0
    f = @(x,Y) crane_ode(x, Y, K0, K1);
    [~, Y] = rk4_system(f, x0, Y0, h, L);
    slope_end = Y(2,end);
end

function dYdx = crane_ode(x, Y, K0, K1)
    % ODE-system för kranprofil Y = [y; y']
    y = Y(1);
    v = Y(2);
    dYdx = [v; -(K0 - K1*x)*y*(1 + v^2)^(3/2)];
end

function [x, Y] = rk4_system(f, x0, Y0, h, L)
    % RK4
    x = x0:h:L;
    N= numel(x);
    Y = zeros(2, N);
    Y(:,1) = Y0;
    for i = 1:N-1
        k1 = f(x(i),       Y(:,i));
        k2 = f(x(i)+h/2,   Y(:,i) + h*k1/2);
        k3 = f(x(i)+h/2,   Y(:,i) + h*k2/2);
        k4 = f(x(i)+h,     Y(:,i) + h*k3);
        Y(:,i+1) = Y(:,i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end

disp(['bestämt K0 = ' num2str(K0_solution, '%.6f')])
disp(['maxhöjd (h) = ' num2str(max_yh, '%.6f')])
disp(['maxhöjd (h/2) = ' num2str(max_yh2, '%.6f')])
disp(['metodfelet= ' num2str(method_err, '%.2e')])
