%% Parametrar och initialvillkor
K0    = 11;
K1    = 0.2;
y0    = 0.1;
s0    = tan(deg2rad(46));  % y'(0)
xspan = [0, 0.5];

%% Toleranser för metodfel
opts_coarse = odeset('RelTol',1e-13,  'AbsTol',1e-14);
opts_fine   = odeset('RelTol',1e-15,  'AbsTol',1e-16);

%% 1) Metodfel via två toleransnivåer
[y_coarse, y_fine, method_err] = method_error(xspan, [y0, s0], K0, K1, opts_coarse, opts_fine);

%% 2) Indatafel via sensitivitetsanalys (±1 %)
input_err = input_error([y0, s0, K0, K1], xspan, opts_coarse);

%% 3) Total felgräns
total_err = method_err + input_err;

%% 4) Visa resultat
disp('Avancerad deluppgift a)')
disp(['  y(0.5)       = ' num2str(y_coarse,'%.8f') ' m'])
disp(['  Metodfel     = ' num2str(method_err,'%.2e') ' m'])
disp(['  Indatafel ≤  = ' num2str(input_err,'%.2e') ' m'])
disp(['  Total fel ≤  = ' num2str(total_err,'%.2e') ' m'])



%% -------- Lokala funktioner --------

function [y1, y2, err] = method_error(xspan, Y0, K0, K1, opts1, opts2)
% Löser ODE med två olika toleranser och returnerar båda lösningarna + fel
    y1 = solve_ode(xspan, Y0, K0, K1, opts1);
    y2 = solve_ode(xspan, Y0, K0, K1, opts2);
    err = abs(y2 - y1);
end

function err = input_error(params, xspan, opts)
% Gör sensitivitetsanalys: ±1 % variation av varje parameter
    base_ref = solve_ode(xspan, params(1:2), params(3), params(4), opts);
    err_contrib = zeros(1,4);
    for j = 1:4
        p = params;
        p(j) = 1.01 * p(j);
        y_p = solve_ode(xspan, p(1:2), p(3), p(4), opts);
        err_contrib(j) = abs(y_p - base_ref);
    end
    err = sum(err_contrib);
end

function y_end = solve_ode(xspan, Y0, K0, K1, opts)
% Kör ode45 för kransystemet och returnerar y vid x = xspan(end)
    odefun = @(x, Y) crane_ode(x, Y, K0, K1);
    [~, Y] = ode45(odefun, xspan, Y0, opts);
    y_end = Y(end,1);
end

function dYdx = crane_ode(x, Y, K0, K1)
% ODE-system för vattenkran: Y = [y; y']
    y = Y(1);
    v = Y(2);
    dYdx = [ v;
            -(K0 - K1*x) * y * (1 + v^2)^(3/2) ];
end
