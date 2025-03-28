function solve_deluppgift2_secant
    % Deluppgift 2: Justera K0 så att y'(0.5) = -0.51 med sekantmetoden
    
    % Parametrar för sekantmetoden
    tol = 1e-6;
    maxIter = 100;
    % Initiala gissningar för K0
    K0_0 = 10.0;
    K0_1 = 11.0;
    
    iter = 0;
    while abs(K0_1 - K0_0) > tol && iter < maxIter
        f0 = f_target(K0_0);
        f1 = f_target(K0_1);
        
        % Sekantformeln
        K0_new = K0_1 - f1 * (K0_1 - K0_0) / (f1 - f0);
        
        % Uppdatera gissningarna
        K0_0 = K0_1;
        K0_1 = K0_new;
        iter = iter + 1;
    end
    
    fprintf('Konvergerat K0: %.12f\n', K0_1);
    fprintf('Antal iterationer: %d\n', iter);
    
    % För att verifiera, integrera med det erhållna K0
    K1 = 0.2;
    y0 = 0.1;
    s0 = tan(deg2rad(46));
    Y0 = [y0; s0];
    x0 = 0;
    L = 0.5;
    h = 0.0001;
    f = @(x, Y) ode_system(x, Y, K0_1, K1);
    [x, Y] = rk4_system(f, x0, Y0, h, L);
    v_end = Y(2, end);
    fprintf('y''(0.5) = %.8f\n', v_end);
end

% Målfunktion: returnerar felet f_target(K0) = y'(0.5) + 0.51
function err = f_target(K0)
    K1 = 0.2;
    y0 = 0.1;
    s0 = tan(deg2rad(46));
    Y0 = [y0; s0];
    x0 = 0;
    L = 0.5;
    h = 0.001;
    f = @(x, Y) ode_system(x, Y, K0, K1);
    [~, Y] = rk4_system(f, x0, Y0, h, L);
    v_end = Y(2, end);
    err = v_end + 0.51;  % Eftersom vi vill att y'(0.5) ska vara -0.51
end

% ODE-systemet
function dYdx = ode_system(x, Y, K0, K1)
    % Y(1) = y, Y(2) = y'
    y = Y(1);
    v = Y(2);
    dYdx = zeros(2,1);
    dYdx(1) = v;
    dYdx(2) = - (K0 - K1 * x) * y * (1 + v^2)^(3/2);
end

% RK4-metod för system av differentialekvationer
function [x, Y] = rk4_system(f, x0, Y0, h, L)
    x = x0:h:L;
    N = length(x);
    Y = zeros(length(Y0), N);
    Y(:,1) = Y0;
    
    for i = 1:N-1
        k1 = f(x(i), Y(:,i));
        k2 = f(x(i) + h/2, Y(:,i) + h * k1/2);
        k3 = f(x(i) + h/2, Y(:,i) + h * k2/2);
        k4 = f(x(i) + h, Y(:,i) + h * k3);
        Y(:,i+1) = Y(:,i) + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
end
