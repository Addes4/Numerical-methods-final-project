% Parametrar
x0           = 0;
L            = 0.5;
y0           = 0.1;
K1           = 0.2;
target_slope = -0.51;
target_height= 0.255;

% ----- Enradiga hjälpfunktioner -----
get_slope = @(K0,s0) ...
    deval( ...
      ode45( @(x,Y)[Y(2); -(K0-K1*x)*Y(1)*(1+Y(2)^2)^(3/2)], ...
             [x0, L], [y0; s0] ), ...
      L, 2 );

get_height = @(K0,s0) ...
    max( ...
      deval( ...
        ode45( @(x,Y)[Y(2); -(K0-K1*x)*Y(1)*(1+Y(2)^2)^(3/2)], ...
               [x0, L], [y0; s0] ), ...
        linspace(x0, L, 500), 1 ) );

% ----- Rotsökning -----
% 1) Hitta startlutning s0 så att maxhöjden blir rätt
find_K0 = @(s0) fzero(@(K0) get_slope(K0,s0) - target_slope, 10);
height_err = @(s0) get_height(find_K0(s0), s0) - target_height;
s0_sol = fzero(height_err, [tan(deg2rad(40)), tan(deg2rad(50))]);

% 2) Hitta K0 med funnet s0
K0_sol = find_K0(s0_sol);

fprintf('Resultat: s0 = %.4f°  ,  K0 = %.6f\n', rad2deg(s0_sol), K0_sol);
