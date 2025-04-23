%parametrar
K0 = 11;
K1 = 0.2;
y0 = 0.1;
s0 = tan(deg2rad(46));
Y0 = [y0; s0];

%steglängder o intervall
x0 = 0;
L  = 0.5;
h  = 1e-4;  % steg
h2 = h/2;   % halvasteg

% använd den generiska ode_system-funktionen istället för crane_ode
f = @(x,Y) ode_system(x, Y, K0, K1);

% kör RK4 för h med extern rk4_system
[~, Yh ]  = rk4_system(f, x0, Y0, h,  L);
yh  = Yh(1,end);

% kör RK4 för h/2
[~, Yh2]  = rk4_system(f, x0, Y0, h2, L);
yh2 = Yh2(1,end);

% metodfel via två steglängder
err_method = abs(yh2 - yh);

disp(['y(0.5) = ' num2str(yh2, '%.6f')])
disp(['metodfel = '    num2str(err_method,  '%.2e')])
