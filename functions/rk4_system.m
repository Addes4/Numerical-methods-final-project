function [x, Y] = rk4_system(f, x0, Y0, h, L)
    % grid
    x = x0:h:L;
    N = numel(x);

    % lösningsmatris
    m = numel(Y0);
    Y = zeros(m, N);
    Y(:,1) = Y0;

    % RK4 iteration
    for i = 1:N-1
        xi = x(i);
        Yi = Y(:, i);
        k1 = f(xi,             Yi);
        k2 = f(xi + h/2,       Yi + (h/2)*k1);
        k3 = f(xi + h/2,       Yi + (h/2)*k2);
        k4 = f(xi + h,         Yi + h*k3);
        Y(:, i+1) = Yi + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end
