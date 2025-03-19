clear; clc; close all;

%% Parametrar
K0 = 11;                 % Konstant K0
K1 = 0.2;                % Konstant K1
K = @(x) K0 - K1*x;      % Värmelednings-/böjningskonstanten K(x)

y0 = 0.1;                          % Starthöjd, y(0)
theta0 = 46 * pi/180;                % Startvinkel i radianer (46°)
yp0 = tan(theta0);                   % Startlutning, y'(0)

%% Differentialekvation
% Vi omvandlar ekvationen y'' = -K(x)*y*(1+(y')^2)^(3/2)
% till ett system av första ordningens ekvationer:
%   Y(1) = y,  Y(2) = y'
%   Y(1)' = y' = Y(2)
%   Y(2)' = y'' = -K(x)*y*(1+(y')^2)^(3/2)
odefun = @(x, Y) [ Y(2); - K(x)*Y(1)*(1 + Y(2)^2)^(3/2) ];

%% Tidsintervall (x-intervallet, där x mäter avståndet i sidled)
xspan = [0, 0.5];

%% Lös ODE med ode45
options = odeset('RelTol',1e-8, 'AbsTol',1e-10);
[x, Y] = ode45(odefun, xspan, [y0; yp0], options);

%% Extrahera lösningen
y_sol = Y(:,1);   % y(x)
yp_sol = Y(:,2);  % y'(x)

%% Sluthöjd vid x = 0.5
final_y = y_sol(end);

%% Plotta kranens profil
figure;
plot(x, y_sol, 'b-', 'LineWidth',2);
xlabel('x [m]');
ylabel('y [m]');
title('Vattenkranens profil (sett från sidan)');
grid on;

%% Skriv ut resultatet
disp(['Vid x = 0.5 slutar kranen på höjden y = ', num2str(final_y, '%.4f')]);
