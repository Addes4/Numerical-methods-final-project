function solve_deluppgift1
    % Parametrar
    K0 = 11;         % Givet värde
    K1 = 0.2;        % Givet värde

    % Initialvillkor
    y0 = 0.1;
    s0 = tan(deg2rad(46));  % y'(0) = tan(46°)
    Y0 = [y0; s0];
    
    % Intervall och steglängd
    x0 = 0;
    L = 0.5;
    h = 0.0001;  % steglängd, välj h litet för noggrannhet

    % Anonym funktion som binder parametrarna K0 och K1
    f = @(x, Y) ode_system(x, Y, K0, K1);
    
    % Lös systemet med RK4
    [x, Y] = rk4_system(f, x0, Y0, h, L);
    
    % Extrahera lösningar
    y = Y(1, end);   % lösning för y vid x = 0.5
    v = Y(2, end);   % lösning för y' vid x = 0.5
    
    % Visa värden vid x = 0.5
    fprintf('y(0.5) = %.8f\n', y);
    fprintf('y''(0.5) = %.8f\n', v);
end

% Funktion som definierar ODE-systemet
function dYdx = ode_system(x, Y, K0, K1)
    % Y(1) = y, Y(2) = y'
    y = Y(1);
    v = Y(2);
    dYdx = zeros(2, 1);
    dYdx(1) = v;
    dYdx(2) = - (K0 - K1 * x) * y * (1 + v^2)^(3/2);
end

% RK4-metod för system av differentialekvationer
function [x, Y] = rk4_system(f, x0, Y0, h, L)
    x = x0:h:L;
    N = length(x);
    Y = zeros(length(Y0), N);
    Y(:, 1) = Y0;
    
    for i = 1:N-1
        k1 = f(x(i), Y(:, i));
        k2 = f(x(i) + h/2, Y(:, i) + h * k1 / 2);
        k3 = f(x(i) + h/2, Y(:, i) + h * k2 / 2);
        k4 = f(x(i) + h, Y(:, i) + h * k3);
        Y(:, i+1) = Y(:, i) + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
end
    
