% Funktion som definierar ODE-systemet
function dYdx = ode_system(x, Y, K0, K1)
    % Y(1) = y, Y(2) = y'
    y = Y(1);
    v = Y(2);
    dYdx = zeros(2, 1);
    dYdx(1) = v;
    dYdx(2) = - (K0 - K1 * x) * y * (1 + v^2)^(3/2);
end
