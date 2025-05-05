% parameter och startvärden
K1 = 0.2;
Y0 = [0.1; tan(deg2rad(46))];   % [y(0); y'(0)]
x0 = 0;
L = 0.5;
h = 1e-5;

% sekantmetod för K0
target_slope = -0.51;
K0_prev = 5;
K0_curr = 15;
f_prev = compute_slope(K0_prev, x0, Y0, h, L, K1) - target_slope;
f_curr = compute_slope(K0_curr, x0, Y0, h, L, K1) - target_slope;

% hjälpfunktioner
function slope_end = compute_slope(K0, x0, Y0, h, L, K1)
    f = @(x,Y) ode_system(x, Y, K0, K1);
    [~, Y] = rk4_system(f, x0, Y0, h, L);
    slope_end = Y(2,end);
end

tol_K = 1e-5;
maxIter = 50;
iter = 0;
while abs(K0_curr - K0_prev) > tol_K && iter < maxIter
    K0_new  = K0_curr - f_curr*(K0_curr - K0_prev)/(f_curr - f_prev);
    K0_prev = K0_curr;  f_prev = f_curr;
    K0_curr = K0_new;
    f_curr  = compute_slope(K0_curr, x0, Y0, h, L, K1) - target_slope;
    iter = iter + 1;
end
K0_solution = K0_curr;

% definiera ODE
f = @(x,Y) ode_system(x, Y, K0_solution, K1);

% RK4 med steglängd h
[x_h,  Yh  ] = rk4_system(f, x0, Y0,  h,  L);
[raw_max_yh, idx_h ] = max(Yh(1,:));

% interpolation för h
if idx_h>1 && idx_h<length(x_h)
    xs_h = x_h(idx_h-1:idx_h+1);
    ys_h = Yh(1, idx_h-1:idx_h+1);
    p_h  = polyfit(xs_h, ys_h, 2);
    x_max_interp_h = -p_h(2)/(2*p_h(1));
    max_yh_interp = polyval(p_h, x_max_interp_h);
else
    x_max_interp_h = x_h(idx_h);
    max_yh_interp  = raw_max_yh;
end

% RK4 med steglängd h/2
h2 = h/2;
[x_h2, Yh2 ] = rk4_system(f, x0, Y0, h2, L);
[raw_max_yh2, idx_h2] = max(Yh2(1,:));

% interpolation för h/2
if idx_h2>1 && idx_h2<length(x_h2)
    xs_h2 = x_h2(idx_h2-1:idx_h2+1);
    ys_h2 = Yh2(1, idx_h2-1:idx_h2+1);
    p_h2 = polyfit(xs_h2, ys_h2, 2);
    x_max_interp_h2 = -p_h2(2)/(2*p_h2(1));
    max_yh2_interp = polyval(p_h2, x_max_interp_h2);
else
    x_max_interp_h2 = x_h2(idx_h2);
    max_yh2_interp  = raw_max_yh2;
end

% metodfel baserat på interpolerade värden
method_err = abs(max_yh2_interp - max_yh_interp);

% utskrift
disp(['K0 = ' num2str(K0_solution, '%.6f')]);
disp(['maxhöjd (h) = ' num2str(raw_max_yh,  '%.6f')]);
disp(['maxhöjd (h/2) = ' num2str(raw_max_yh2, '%.6f')]);
disp(['metodfel (interpolerat) = ' num2str(method_err,  '%.2e')]);
