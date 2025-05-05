% parametrar
K0 = 10.675321;  % från uppgift b
K1 = 0.2;
x0 = 0; L = 0.5;
Y0 = [0.1; tan(deg2rad(46))];
h  = 1e-5;

% RK4 med h
f = @(x,Y) ode_system(x,Y,K0,K1);
[x1,Y1] = rk4_system(f, x0, Y0, h,  L);

% RK4 med h/2
[x2,Y2] = rk4_system(f, x0, Y0, h/2, L);


function [x_max,y_max] = quadratic_interp(x,y,idx)
    n = numel(x);
    if idx>1 && idx<n
        xm = x(idx-1:idx+1).';    % kolumnvektor
        ym = y(idx-1:idx+1).';    % kolumnvektor
        A  = [xm.^2, xm, ones(3,1)];
        c  = A \ ym;              % [a; b; c]
        x_max = -c(2)/(2*c(1));   % vertex
        y_max = c(1)*x_max^2 + c(2)*x_max + c(3);
    else
        x_max = x(idx);
        y_max = y(idx);
    end
end

% hitta och interpolera max
[~,i1] = max(Y1(1,:));
[~,y1] = quadratic_interp(x1, Y1(1,:), i1);
[~,i2] = max(Y2(1,:));
[~,y2] = quadratic_interp(x2, Y2(1,:), i2);

% metodfel
err = abs(y2 - y1);

fprintf('Maxhöjd   = %.6f\n', y1);
fprintf('Metodfel  = %.2e\n', err);
