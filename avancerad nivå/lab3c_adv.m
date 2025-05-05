% parametrar
K0 = 10.675321;
K1 = 0.2;
Y0 = [0.1; tan(deg2rad(46))];
xspan = [0, 0.5];

% tol för ode45
opts_work = odeset('RelTol',1e-12,'AbsTol',1e-13);
opts_ref  = odeset('RelTol',1e-14,'AbsTol',1e-15);

% ODE-funktion
odeFun = @(x,y)[y(2);-(K0 - K1*x)*y(1)*(1 + y(2)^2)^(3/2)];

% arbets- och referenstolerans
[x_w, Y_w] = ode45(odeFun, xspan, Y0, opts_work);
[x_r, Y_r] = ode45(odeFun, xspan, Y0, opts_ref);

%kvad interp
function [ymax, xmax] = interp_quad(x, y)
    [~, idx] = max(y);
    if idx>1 && idx<numel(x)
        xs = x(idx-1:idx+1);
        ys = y(idx-1:idx+1);
        p  = polyfit(xs, ys, 2);
        xmax = -p(2)/(2*p(1));
        ymax = polyval(p, xmax);
    else
        xmax = x(idx);
        ymax = y(idx);
    end
end

% hitta maxhöjd m kvad interp
[ max_w, xw ] = interp_quad(x_w, Y_w(:,1));
[ max_r, xr ] = interp_quad(x_r, Y_r(:,1));

% metodfelet
method_err = abs(max_r - max_w);

fprintf('K0   = %.6f\n', K0);
fprintf('maxhöjd  = %.6f at x = %.6f\n', max_w, xw);
fprintf('metodfelet  = %.2e\n', method_err);
