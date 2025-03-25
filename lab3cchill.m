clear; clc; close all;

%% Parametrar för vattenkranproblemet
L = 0.5;                % integrationsintervall: x i [0, L]
K0 = 10.8;              % justerat värde (från deluppgift 2)
K1 = 0.2;
y0 = 0.1;
s0 = tan(deg2rad(46));   % y'(0) = tan(46°)
Y0 = [y0; s0];

%% Första beräkningen med n1 steg
n1 = 500;
h1 = L/n1;
x1 = linspace(0, L, n1+1);
[~, Y1] = rk4_system(@(x,Y) odefun(x,Y,K0,K1), x1, Y0);
% Y1(1,:) = y; Y1(2,:) = y'
index1 = find_sign_change(Y1(2,:));
[xmax1, ymax1] = interpolate_max(x1, Y1(1,:), Y1(2,:), index1);

%% Andra beräkningen med n2 steg
n2 = 1000;
h2 = L/n2;
x2 = linspace(0, L, n2+1);
[~, Y2] = rk4_system(@(x,Y) odefun(x,Y,K0,K1), x2, Y0);
index2 = find_sign_change(Y2(2,:));
[xmax2, ymax2] = interpolate_max(x2, Y2(1,:), Y2(2,:), index2);

%% Richardson-extrapolering
y_extrap = ymax2 + (ymax2 - ymax1)/3;
err_est = abs(y_extrap - ymax2);

%% Utskrift och plot
fprintf('Med n1 = %d: y_max = %.4f\n', n1, ymax1);
fprintf('Med n2 = %d: y_max = %.4f\n', n2, ymax2);
fprintf('Richardson extrapolerat y_max = %.4f\n', y_extrap);
fprintf('Uppskattat fel = %.4e\n', err_est);

figure;
plot(x2, Y2(1,:), 'b-o', 'MarkerSize',4, 'LineWidth',1.5);
hold on;
plot(xmax2, ymax2, 'rx', 'MarkerSize',12, 'LineWidth',2);
xlabel('x'); ylabel('y');
title('Vattenkranens form - maximal höjd');
grid on;

%% Funktioner

function dYdx = odefun(x, Y, K0, K1)
    % ODE-systemet:
    % Y(1) = y, Y(2) = y'
    y = Y(1);
    v = Y(2);
    dYdx = zeros(2,1);
    dYdx(1) = v;
    dYdx(2) = - (K0 - K1*x) * y * (1 + v^2)^(3/2);
end

function [x, Y] = rk4_system(f, x, Y0)
    % RK4-metod för system av differentialekvationer.
    % x: vektor med x-värden, Y0: initialtillstånd (kolonnvektor)
    N = length(x);
    Y = zeros(length(Y0), N);
    Y(:,1) = Y0;
    h = x(2) - x(1);
    for i = 1:N-1
        k1 = f(x(i), Y(:,i));
        k2 = f(x(i) + h/2, Y(:,i) + h*k1/2);
        k3 = f(x(i) + h/2, Y(:,i) + h*k2/2);
        k4 = f(x(i) + h, Y(:,i) + h*k3);
        Y(:,i+1) = Y(:,i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end

function idx = find_sign_change(v)
    % Hitta första index där v byter från positivt till icke-positivt (t.ex. från >0 till <=0)
    idx = [];
    for i = 1:length(v)-1
        if v(i) > 0 && v(i+1) <= 0clear; clc; close all;

%% Parametrar för vattenkranproblemet
L = 0.5;                % integrationsintervall: x i [0, L]
K0 = 10.8;              % justerat värde (från deluppgift 2)
K1 = 0.2;
y0 = 0.1;
s0 = tan(deg2rad(46));   % y'(0) = tan(46°)
Y0 = [y0; s0];

%% Första beräkningen med n1 steg
n1 = 500;
h1 = L/n1;
x1 = linspace(0, L, n1+1);
[~, Y1] = rk4_system(@(x,Y) odefun(x,Y,K0,K1), x1, Y0);
% Y1(1,:) = y; Y1(2,:) = y'
index1 = find_sign_change(Y1(2,:));
[xmax1, ymax1] = interpolate_max(x1, Y1(1,:), Y1(2,:), index1);

%% Andra beräkningen med n2 steg
n2 = 1000;
h2 = L/n2;
x2 = linspace(0, L, n2+1);
[~, Y2] = rk4_system(@(x,Y) odefun(x,Y,K0,K1), x2, Y0);
index2 = find_sign_change(Y2(2,:));
[xmax2, ymax2] = interpolate_max(x2, Y2(1,:), Y2(2,:), index2);

%% Richardson-extrapolering
y_extrap = ymax2 + (ymax2 - ymax1)/3;
err_est = abs(y_extrap - ymax2);

%% Utskrift och plot
fprintf('Med n1 = %d: y_max = %.4f\n', n1, ymax1);
fprintf('Med n2 = %d: y_max = %.4f\n', n2, ymax2);
fprintf('Richardson extrapolerat y_max = %.4f\n', y_extrap);
fprintf('Uppskattat fel = %.4e\n', err_est);

figure;
plot(x2, Y2(1,:), 'b-o', 'MarkerSize',4, 'LineWidth',1.5);
hold on;
plot(xmax2, ymax2, 'rx', 'MarkerSize',12, 'LineWidth',2);
xlabel('x'); ylabel('y');
title('Vattenkranens form - maximal höjd');
grid on;

%% Funktioner

function dYdx = odefun(x, Y, K0, K1)
    % ODE-systemet:
    % Y(1) = y, Y(2) = y'
    y = Y(1);
    v = Y(2);
    dYdx = zeros(2,1);
    dYdx(1) = v;
    dYdx(2) = - (K0 - K1*x) * y * (1 + v^2)^(3/2);
end

function [x, Y] = rk4_system(f, x, Y0)
    % RK4-metod för system av differentialekvationer.
    % x: vektor med x-värden, Y0: initialtillstånd (kolonnvektor)
    N = length(x);
    Y = zeros(length(Y0), N);
    Y(:,1) = Y0;
    h = x(2) - x(1);
    for i = 1:N-1
        k1 = f(x(i), Y(:,i));
        k2 = f(x(i) + h/2, Y(:,i) + h*k1/2);
        k3 = f(x(i) + h/2, Y(:,i) + h*k2/2);
        k4 = f(x(i) + h, Y(:,i) + h*k3);
        Y(:,i+1) = Y(:,i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end

function idx = find_sign_change(v)
    % Hitta första index där v byter från positivt till icke-positivt (t.ex. från >0 till <=0)
    idx = [];
    for i = 1:length(v)-1
        if v(i) > 0 && v(i+1) <= 0
            idx = i;
            break;
        end
    end
    if isempty(idx)
        error('Ingen signalförändring (nollgenomgång) hittades i lutningen.');
    end
end

function [x_max, y_max] = interpolate_max(x, y, dydx, idx)
    % Linjär interpolation mellan x(idx) och x(idx+1) för att finna x där y' = 0.
    % Förutsätter att dydx(idx) > 0 och dydx(idx+1) < 0.
    x1 = x(idx);
    x2 = x(idx+1);
    v1 = dydx(idx);
    v2 = dydx(idx+1);
    % Interpolera för x_max:
    x_max = x1 - v1*(x2 - x1)/(v2 - v1);
    % Interpolera för y_max:
    y_max = y(idx) + (y(idx+1) - y(idx))*(x_max - x1)/(x2 - x1);
end
            idx = i;
            break;
        end
    end
    if isempty(idx)
        error('Ingen signalförändring (nollgenomgång) hittades i lutningen.');
    end
end

function [x_max, y_max] = interpolate_max(x, y, dydx, idx)
    % Linjär interpolation mellan x(idx) och x(idx+1) för att finna x där y' = 0.
    % Förutsätter att dydx(idx) > 0 och dydx(idx+1) < 0.
    x1 = x(idx);
    x2 = x(idx+1);
    v1 = dydx(idx);
    v2 = dydx(idx+1);
    % Interpolera för x_max:
    x_max = x1 - v1*(x2 - x1)/(v2 - v1);
    % Interpolera för y_max:
    y_max = y(idx) + (y(idx+1) - y(idx))*(x_max - x1)/(x2 - x1);
end
