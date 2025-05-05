% parameter 0 startvillkor
K1 = 0.2;
Y0 = [0.1; tan(deg2rad(46))];  % [y(0); y'(0)]
xspan = [0, 0.5];

% sluttlutning y'(L)
target_slope = -0.51;

% inställ ode45
opts_work = odeset('RelTol',1e-12,'AbsTol',1e-13);
opts_ref  = odeset('RelTol',1e-14,'AbsTol',1e-15);

% K0 via fzero + inline ode45/deval ---
g = @(K0) deval(ode45(@(x,y)[y(2);-(K0 - K1*x)*y(1)*(1 + y(2)^2)^(3/2)], ...
           xspan, Y0, opts_work),xspan(end), 2) - target_slope;

K0_initial  = 11;
K0_solution = fzero(g, K0_initial);

% def ode-system m funktionshandle
odeFun = @(x,y)[y(2);-(K0_solution - K1*x)*y(1)*(1 + y(2)^2)^(3/2)];

%integrera m båda inst
[x_w, Y_w] = ode45(odeFun, xspan, Y0, opts_work);
[x_r, Y_r] = ode45(odeFun, xspan, Y0, opts_ref);

%rå maximihöjd 
[raw_max_yw, idx_w] = max(Y_w(:,1));
[raw_max_yr, idx_r] = max(Y_r(:,1));

% 3e ordnings interp
if idx_w > 2 && idx_w < length(x_w)-1
    x_seg = x_w(idx_w-2:idx_w+1);
    y_seg = Y_w(idx_w-2:idx_w+1,1);
    p = polyfit(x_seg, y_seg, 3);  % 3e ordnings polynom
    dp = polyder(p);
    r_all = roots(dp);  % Kritiska punkter
    % reella rötter inom segmentet
    r = r_all(imag(r_all)==0 & r_all>=min(x_seg) & r_all<=max(x_seg));
    if ~isempty(r)
        %den som ger högst polynomvärde
        [~, imax]      = max(polyval(p, r));
        x_max_interp   = r(imax);
        max_yw_interp  = polyval(p, x_max_interp);
    else
        x_max_interp  = x_w(idx_w);
        max_yw_interp = raw_max_yw;
    end
else
    % kantfall använd råvärdet
    x_max_interp  = x_w(idx_w);
    max_yw_interp = raw_max_yw;
end

%fel
method_err = abs(raw_max_yr - raw_max_yw);

fprintf('K0 = %.6f\n', K0_solution);
fprintf('rå maxhöjd = %.6f\n', raw_max_yw);
fprintf('interp maxhöjd (3a) = %.6f\n', max_yw_interp);
fprintf('maxhöjd (h/2) (rå)= %.6f\n', raw_max_yr);
fprintf('metodfelet= %.2e\n', method_err);
