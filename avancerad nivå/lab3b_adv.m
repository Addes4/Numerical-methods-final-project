
%parametrar
K1 = 0.2;
Y0 = [0.1; tan(deg2rad(46))];  % [y(0); y'(0)]
xspan = [0, 0.5];  % x-intervall

target_slope = -0.51;  % mål slutlutning

%solverinställningar för metodfel
opts_work = odeset('RelTol',1e-12, 'AbsTol',1e-13);
opts_ref = odeset('RelTol',1e-14, 'AbsTol',1e-15);

%  funktionshandle för att beräkna slutlutning
slope_work = @(K0) compute_slope(K0, xspan, Y0, K1, opts_work);
slope_ref  = @(K0) compute_slope(K0, xspan, Y0, K1, opts_ref);

%rot-funktion för fzero (baserat på arbetslösningen)
g = @(K0) slope_work(K0) - target_slope;

%hitta K0 med fzero
K0_initial = 11;
K0_solution = fzero(g, K0_initial);

%beräkna metodfel vid hittat K0
yw = slope_work(K0_solution);
yr = slope_ref(K0_solution);
method_err = abs(yr - yw);

function slope_end = compute_slope(K0, xspan, Y0, K1, opts)
    %kör ODE och returnerar y'(xspan(end))
    odeFun = @(x, Y) [Y(2);-(K0 - K1*x)*Y(1)*(1 + Y(2)^2)^(3/2)];
    [~, Y] = ode45(odeFun, xspan, Y0, opts);
    slope_end = Y(end,2);
end

disp(['nya K0 = ' num2str(K0_solution,'%.8f')])
disp(['y''(0.5) = ' num2str(yw,'%.6f')])
disp(['metodfelet= ' num2str(method_err,'%.2e')])
