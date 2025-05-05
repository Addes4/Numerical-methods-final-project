
% parameter och startvärd
K1 = 0.2;
Y0 = [0.1; tan(deg2rad(46))];   % [y(0); y'(0)]
x0 = 0;
L = 0.5;
h = 1e-5;

%  K0 med sekantmetod
target_slope = -0.51;
K0_prev = 5;  % gissning 1
K0_curr = 15; % gissning 2
f_prev = compute_slope(K0_prev, x0, Y0, h, L, K1) - target_slope;
f_curr = compute_slope(K0_curr, x0, Y0, h, L, K1) - target_slope;

tol_K = 1e-5;
maxIter = 50;
iter = 0;
while abs(K0_curr - K0_prev) > tol_K && iter < maxIter
    K0_new  = K0_curr - f_curr * (K0_curr - K0_prev) / (f_curr - f_prev);
    K0_prev = K0_curr;  f_prev = f_curr;
    K0_curr = K0_new;
    f_curr  = compute_slope(K0_curr, x0, Y0, h, L, K1) - target_slope;
    iter = iter + 1;
end
K0_solution = K0_curr; % lösning K0

%  kranprofil och maximihöjd
f = @(x,Y) ode_system(x, Y, K0_solution, K1);
% körning m steglängd h
[x_h, Yh ] = rk4_system(f, x0, Y0, h,    L);
% grov maximihöjd o index
[raw_max_yh, idx] = max(Yh(1,:));
x_vals = x_h;  % X-vektor fr RK4

% 2a grads interpolering
if idx > 1 && idx < length(x_vals)
    x_segment = x_vals(idx-1:idx+1);
    y_segment = Yh(1, idx-1:idx+1);
    p = polyfit(x_segment, y_segment, 2); % 2 ordningens polynom
    x_max_interp = -p(2) / (2*p(1));
    max_yh_interp = polyval(p, x_max_interp);
else
    % om maximalpunkt är i kanten, använd råvärde
    x_max_interp = x_vals(idx);
    max_yh_interp = raw_max_yh;
end

%  compute_slope 
function slope_end = compute_slope(K0, x0, Y0, h, L, K1)
    f = @(x,Y) ode_system(x, Y, K0, K1);
    [~, Y] = rk4_system(f, x0, Y0, h, L);
    slope_end = Y(2,end);
end

%  steglängd h/2 för feluppskattning
h2 = h/2;
[~, Yh2 ] = rk4_system(f, x0, Y0, h2,   L);
max_yh2 = max(Yh2(1,:));

method_err = abs(max_yh2 - max_yh);

disp([' K0 = ' num2str(K0_solution, '%.6f')]);
disp(['maxhöjd (h) = ' num2str(raw_max_yh,     '%.6f')]);
disp(['interpolerad maxhöjd = ' num2str(max_yh_interp, '%.6f')]);
disp(['maxhöjd (h/2) = ' num2str(max_yh2,        '%.6f')]);
disp(['metodfel (råa) = ' num2str(method_err,     '%.2e')]);
