
%parametrar
K0 = 11;
K1 = 0.2;
y0= 0.1;
s0 = tan(deg2rad(46));
Y0 = [y0; s0];

%stegländer o interv
x0 = 0;
L = 0.5;
h = 1e-4; %steg
h2= h/2; %halvasteg

%kör RK4 för h
f = @(x,Y) crane_ode(x,Y,K0,K1);
[~, Yh ] = rk4_system(f, x0, Y0, h,  L);
yh = Yh(1,end);

%kör RK4 för h/2
[~, Yh2] = rk4_system(f, x0, Y0, h2, L);
yh2 = Yh2(1,end);

% felet
err_method = abs(yh2 - yh);

% funktioner
function dYdx = crane_ode(x, Y, K0, K1)
    % ODE-system för vattenkran: Y = [y; y']
    y = Y(1);
    v = Y(2);
    dYdx = [ v;
    -(K0 - K1*x) * y * (1 + v^2)^(3/2) ];
end

function [x, Y] = rk4_system(f, x0, Y0, h, L)
    %RK4
    x = x0:h:L;
    N = numel(x);
    Y = zeros(length(Y0), N);
    Y(:,1) = Y0;
    for i = 1:N-1
        k1 = f(x(i),  Y(:,i));
        k2 = f(x(i)+h/2, Y(:,i) + h*k1/2);
        k3 = f(x(i)+h/2, Y(:,i) + h*k2/2);
        k4 = f(x(i)+h, Y(:,i) + h*k3);
        Y(:,i+1) = Y(:,i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end

disp(['y(0.5) = ' num2str(yh2,'%.6f')])
disp(['metodfel= ' num2str(err_method,'%.2e')])
