
%parametrar
K0= 11;
K1 = 0.2;
Y0 = [0.1; tan(deg2rad(46))]; %startvillkor [y(0); y'(0)]
xspan = [0, 0.5]; % x-intervall

% solverinställng för metodfel
opts_work = odeset('RelTol',1e-12,  'AbsTol',1e-13);
opts_ref = odeset('RelTol',1e-14,  'AbsTol',1e-15);

%beräkna arbets och referenslösning
y_work = solve_crane(xspan, Y0, K0, K1, opts_work);
y_ref = solve_crane(xspan, Y0, K0, K1, opts_ref);

%felet
method_err = abs(y_ref - y_work);

%  funktioner
function y_end = solve_crane(xspan, Y0, K0, K1, opts)
    %kör ODE med ode45 och returnera y då x = xspan(end)
    odeFun = @(x, Y) [Y(2);
    -(K0 - K1*x)*Y(1)*(1 + Y(2)^2)^(3/2)];
    [~, Y] = ode45(odeFun, xspan, Y0, opts);
    y_end = Y(end,1);
end

disp(['y(0.5)= ' num2str(y_work,'%.8f')])
disp(['metodfel= ' num2str(method_err,'%.2e')])
