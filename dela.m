clear; clc; close all;

%% Parametrar
L = 3.60;    % Stavens längd [m]
T0 = 310;    % Temperatur vid x = 0 [K]
TL = 450;    % Temperatur vid x = L [K]

%% Nät 1: n1 = 576
n1 = 576;
h1 = L/n1;
x1 = linspace(0, L, n1+1)';  % Kolumnvektor med n1+1 noder
index1 = round(1.65/h1) + 1; % Index för x = 1.65
T_1 = rvp(L, T0, TL, x1, h1);
T_1_165 = T_1(index1);

%% Nät 2: n2 = 1152
n2 = 1152;
h2 = L/n2;
x2 = linspace(0, L, n2+1)';
index2 = round(1.65/h2) + 1;
T_2 = rvp(L, T0, TL, x2, h2);
T_2_165 = T_2(index2);

%% Richardson-extrapolation
T_extrap = T_2_165 + (T_2_165 - T_1_165)/3;
err_est = abs(T_extrap - T_2_165);

%% Utskrift
disp(['T(1.65) med n1 = ', num2str(T_1_165, '%.4f')]);
disp(['T(1.65) med n2 = ', num2str(T_2_165, '%.4f')]);
disp(['Richardson-extr. T(1.65) = ', num2str(T_extrap, '%.4f')]);
disp(['Uppskattat fel = ', num2str(err_est, '%.4e')]);

%% Plottning
figure;
plot(x2, T_2, 'b-o');
xlabel('x [m]');
ylabel('T(x) [K]');
title('Temperaturfördelning i staven');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Funktion: rvp
function T = rvp(L, T0, TL, x, h)
    % rvp löser randvärdesproblemet med finite differensmetoden.
    % Indata:
    %   L   - Stavens längd
    %   T0  - Temperatur vid x = 0
    %   TL  - Temperatur vid x = L
    %   x   - Kolumnvektor med noder i intervallet [0, L]
    %   h   - Steglängd (xₖ₊₁ − xₖ)
    %
    % Utdata:
    %   T   - Kolumnvektor med approximativa temperaturvärden vid noderna.
    
    N = length(x) - 1;       % Totalt antal intervall
    N_int = N - 1;           % Antal inre noder (utan randvärden)
    
    % Initiera matris och högerled
    A = zeros(N_int, N_int);
    b = zeros(N_int, 1);
    
    % Loopa över de inre noderna
    for i = 1:N_int
        idx = i + 1;     % Aktuell nod i den fullständiga x-vektorn
        xi = x(idx);
        
        % Beräkna värdet av k(x) vid halva stegen
        k_ip = 3 + ((x(idx) + x(idx+1))/2)/7;  % k(x_{i+1/2})
        k_im = 3 + ((x(idx-1) + x(idx))/2)/7;    % k(x_{i-1/2})
        
        % Sätt upp koefficienterna för noden idx
        A(i,i) = (k_ip + k_im) / h^2;
        if i > 1
            A(i,i-1) = - k_im / h^2;
        else
            % Justera b med T0 för den första inre noden
            b(i) = b(i) + k_im * T0 / h^2;
        end
        if i < N_int
            A(i,i+1) = - k_ip / h^2;
        else
            % Justera b med TL för den sista inre noden
            b(i) = b(i) + k_ip * TL / h^2;
        end
        
        % Lägg till effekten av värmekällan Q(x) = 280*exp(- (x - L/2)^2)
        b(i) = b(i) + 280 * exp(- (xi - L/2)^2);
    end
    
    % Lös det linjära systemet A*T_interior = b
    T_interior = A\b;
    T = [T0; T_interior; TL];
end
