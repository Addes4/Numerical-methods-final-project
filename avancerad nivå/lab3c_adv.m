% parameter 0 startvillkor
K1 = 0.2;
Y0 = [0.1; tan(deg2rad(46))];
xspan = [0, 0.5];
target_slope = -0.51;

% inställningar för ode45
opts_work = odeset('RelTol',1e-12,'AbsTol',1e-13);
opts_ref  = odeset('RelTol',1e-14,'AbsTol',1e-15);

% hitta K0 med fzero
g = @(K0) deval(ode45(@(x,y)[y(2); -(K0 - K1*x)*y(1)*(1 + y(2)^2)^(3/2)], ...
           xspan, Y0, opts_work), xspan(end), 2) - target_slope;
K0_solution = fzero(g, 11);

% definiera ODE-system
odeFun = @(x,y)[y(2); -(K0_solution - K1*x)*y(1)*(1 + y(2)^2)^(3/2)];

% lös m ode45 och 2 toleranser
[x_w, Y_w] = ode45(odeFun, xspan, Y0, opts_work);
[x_r, Y_r] = ode45(odeFun, xspan, Y0, opts_ref);

% 3e ordningen interp arbetslösning
[max_yw_interp, ~] = interp_poly3(x_w, Y_w(:,1));

% 3e ordningen interp  referenslösning
[max_yr_interp, ~] = interp_poly3(x_r, Y_r(:,1));

% metodfel
method_err = abs(max_yr_interp - max_yw_interp);

fprintf('K0 = %.6f\n', K0_solution);
fprintf('maxhöjd (interp work) = %.6f\n', max_yw_interp);
fprintf('maxhöjd (interp ref)  = %.6f\n', max_yr_interp);
fprintf('metodfelet = %.2e\n', method_err);

%interp funktion
function [ymax, xmax] = interp_poly3(x, y)
    [~, idx] = max(y);
    if idx > 2 && idx < length(x)-1
        x_seg = x(idx-2:idx+1);
        y_seg = y(idx-2:idx+1);
        p = polyfit(x_seg, y_seg, 3);
        dp = polyder(p);
        r = roots(dp);
        r = r(imag(r)==0 & r>=min(x_seg) & r<=max(x_seg));
        if ~isempty(r)
            [ymax, i] = max(polyval(p, r));
            xmax = r(i);
            return;
        end
    end
    xmax = x(idx);
    ymax = y(idx);
end
