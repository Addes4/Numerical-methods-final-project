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