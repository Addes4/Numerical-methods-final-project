
% parametrar
K0 = 11;
K1 = 0.2;
y0 = 0.1;
slope0= tan(deg2rad(46));  % startlutning y'(0)

% lös intervall
x0 = 0;
L = 0.5;

% Initieringsvektor
Y0 = [y0; slope0];

% antal steg för huvudkörning
N = 100000;

% Kör RK4 med steglängd h
Y_h = rk4(@craneODE, x0, L, Y0, N,   K0, K1);
y_end_h = Y_h(1,end);

% Kör RK4 med steglängd h/2 för felskattning
Y_h2 = rk4(@craneODE, x0, L, Y0, 2*N, K0, K1);
y_end_h2 = Y_h2(1,end);

% Felskattning numerisk metod (global fel = ca skillnad)
err_method = abs(y_end_h2 - y_end_h);

% Osäkerhet i indata +- 1% -> sensitivitet
params_nom = [y0, slope0, K0, K1];
err_inputs = zeros(1,4);
for j = 1:4
    p = params_nom;
    p(j) = p(j)*1.01;             % 1% ökning
    Y0p   = [p(1); p(2)];
    Yp    = rk4(@craneODE, x0, L, Y0p, N, p(3), p(4));
    err_inputs(j) = abs(Yp(1,end) - y_end_h);
end
err_input_bound = sum(err_inputs);

% tot felgräns
err_total = err_method + err_input_bound;

function dY = craneODE(x, Y, K0, K1)
    y  = Y(1);
    dy = Y(2);
    K  = K0 - K1*x;
    d2y = -K * y * (1 + dy^2)^(3/2);
    dY = [dy; d2y];
end

function Y = rk4(odefun, x0, xf, Y0, N, K0, K1)
    h = (xf - x0)/N;
    Y = zeros(2, N+1);
    x = x0;
    Y(:,1) = Y0;
    for i = 1:N
        k1 = odefun(x,         Y(:,i),       K0, K1);
        k2 = odefun(x + h/2,   Y(:,i) + h*k1/2, K0, K1);
        k3 = odefun(x + h/2,   Y(:,i) + h*k2/2, K0, K1);
        k4 = odefun(x + h,     Y(:,i) + h*k3,   K0, K1);
        Y(:,i+1) = Y(:,i) + h*(k1 + 2*k2 + 2*k3 + k4)/6;
        x = x + h;
    end
end

% print
disp('Subuppgift a:')
disp(['  y(' num2str(L, '%.2f') ') = ' num2str(y_end_h2, '%.6f') ' m'])
disp(['  Metodfel ≈ '    num2str(err_method,    '%.2e')    ' m'])
disp(['  Indatafel ≤ '   num2str(err_input_bound,'%.2e')    ' m'])
disp(['  Total fel ≤ '   num2str(err_total,      '%.2e')    ' m'])
