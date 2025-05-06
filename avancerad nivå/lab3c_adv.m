% parametrar
K1 = 0.2;
K0 = 10.6753209;
y0 = 0.1;
s0 = tan(deg2rad(46));
Y0 = [y0; s0];
x0 = 0;
L = 0.5;
h = 1e-3;

% ode funktion
odeFunc = @(x, Y) [Y(2); -(K0 - K1*x)*Y(1)*(1 + Y(2)^2)^(3/2)];

% lös ode system
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[x_vals, Y] = ode45(odeFunc, x0:h:L, Y0, opts);
y_vals = Y(:, 1);

% hitta grovt max
[~, idx] = max(y_vals);

% punkter runt max
x1 = x_vals(idx-2); y1 = y_vals(idx-2);
x2 = x_vals(idx-1); y2 = y_vals(idx-1);
x3 = x_vals(idx);   y3 = y_vals(idx);
x4 = x_vals(idx+1); y4 = y_vals(idx+1);

% 2a grads pol
p_quad = polyfit([x2, x3, x4], [y2, y3, y4], 2);
x_max_quad = -p_quad(2)/(2*p_quad(1));
y_max_quad = polyval(p_quad, x_max_quad);

% 3e grads pol
p_cubic = polyfit([x1, x2, x3, x4], [y1, y2, y3, y4], 3);
dp_cubic = polyder(p_cubic);                      

%fel
error_y = abs(y_max_cubic - y_max_quad);

fprintf('Maxhöjd     = %.8f\n', y_max_quad);
fprintf('Felskattning = %.3e\n', error_y);
