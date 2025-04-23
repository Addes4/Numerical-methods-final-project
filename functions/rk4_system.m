function [x, Y] = rk4_system(f, x0, Y0, h, L)
% RK4_SYSTEM Generic fourth-order Runge-Kutta integrator for ODE systems
%   [x, Y] = rk4_system(f, x0, Y0, h, x_end) integrates Y' = f(x,Y) from x0 to x_end
%   using step size h.
%
% INPUT:
%   f      - function handle @(x, Y) returning column vector dY/dx
%   x0     - initial value of independent variable
%   Y0     - initial column vector of dependent variables
%   h      - integration step size
%   x_end  - final value of independent variable
%
% OUTPUT:
%   x      - row vector of grid points from x0 to x_end with spacing h
%   Y      - matrix of size (length(Y0)×length(x)), each column Y(:,i) is solution at x(i)

    % Create grid
    x = x0:h:L;
    N = numel(x);

    % Initialize solution array
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
