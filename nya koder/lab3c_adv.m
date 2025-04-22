
% parametrar
K1= 0.2;
Y0 = [0.1; tan(deg2rad(46))];  % [y(0); y'(0)]
xspan = [0, 0.5];

%finder K0 från uppgift b med samma inställningar
target_slope = -0.51;
opts_work = odeset('RelTol',1e-12,  'AbsTol',1e-13);
opts_ref = odeset('RelTol',1e-14,  'AbsTol',1e-15);

slope_work = @(K0) compute_slope(K0, xspan, Y0, K1, opts_work);
g = @(K0) slope_work(K0) - target_slope;
K0_initial = 11;
K0_solution= fzero(g, K0_initial);

%   beräkna arbets och referenslösning
[x_w, Y_w] = ode45(@(x,Y) crane_ode(x,Y,K0_solution,K1), xspan, Y0, opts_work);
[x_r, Y_r] = ode45(@(x,Y) crane_ode(x,Y,K0_solution,K1), xspan, Y0, opts_ref);

%hitta maximal höjd
[max_yw, idx_w] = max(Y_w(:,1));    % arbetsprofil
[max_yr, idx_r] = max(Y_r(:,1));    % referensprofil

%felet
method_err = abs(max_yr - max_yw);

function slope_end = compute_slope(K0, xspan, Y0, K1, opts)
    odeFun = @(x,Y) crane_ode(x, Y, K0, K1);
    [~, Y] = ode45(odeFun, xspan, Y0, opts);
    slope_end = Y(end,2);
end

function dYdx = crane_ode(x, Y, K0, K1)
    y = Y(1);
    v = Y(2);
    dYdx = [v; -(K0 - K1*x)*y*(1 + v^2)^(3/2)];
end

disp(['maxhöjd = ' num2str(max_yw, '%.8f')])
disp(['metodfel = ' num2str(method_err, '%.2e')])
