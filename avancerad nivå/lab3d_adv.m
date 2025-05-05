% parametrar
x0 = 0;
L = 0.5;
y0 = 0.1;
K1 = 0.2;
target_slope = -0.51;
target_height= 0.255;

% hjälpfunktion beräkning sluttlutning y'(L)
get_slope = @(K0,s0)deval(ode45(@(x,Y)[Y(2); ...
    -(K0-K1*x)*Y(1)*(1+Y(2)^2)^(3/2)],[x0, L], [y0; s0] ), L, 2 );

% hjälpfunktion beräkning maxhöjd y_max
get_height = @(K0,s0)max(deval(ode45( @(x,Y)[Y(2); ...
    -(K0-K1*x)*Y(1)*(1+Y(2)^2)^(3/2)],[x0, L], [y0; s0] ), ...
        linspace(x0, L, 500), 1 ) );

% hitta s0 så att maxhöjd blir rätt
find_K0 = @(s0)fzero(@(K0) get_slope(K0,s0) - target_slope, 10);
height_err = @(s0)get_height(find_K0(s0), s0) - target_height;
s0_sol = fzero(height_err, [tan(deg2rad(40)), tan(deg2rad(50))]);

% hitta K0 med s0
K0_sol = find_K0(s0_sol); % K0 som ger rätt sluttlutning

fprintf('Resultat: s0 = %.4f°  ,  K0 = %.6f\n', rad2deg(s0_sol), K0_sol);

% 2 toleranser för ode45
opts1 = odeset('RelTol',1e-11,  'AbsTol',1e-12);
opts2 = odeset('RelTol',1e-13, 'AbsTol',1e-14);

% spara ODE-funktionen
odeFun = @(x,Y,K0) [Y(2); -(K0-K1*x)*Y(1)*(1+Y(2)^2)^(3/2)];

% m 1a toleransen
sol1 = ode45(@(x,Y) odeFun(x,Y,K0_sol), [x0, L], [y0; tan(s0_sol)], opts1);
s1   = deval(sol1, L, 2);
Yv1  = deval(sol1, linspace(x0, L, 500));
h1   = max(Yv1(1,:));

% m 2a toleransen
sol2 = ode45(@(x,Y) odeFun(x,Y,K0_sol), [x0, L], [y0; tan(s0_sol)], opts2);
s2   = deval(sol2, L, 2);
Yv2  = deval(sol2, linspace(x0, L, 500));
h2   = max(Yv2(1,:));

% metodfelet
err_slope  = abs(s2 - s1);
err_height = abs(h2 - h1);

fprintf('Metodfel: err_slope = %.2e, err_height = %.2e\n', err_slope, err_height);
