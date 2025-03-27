clear; clc; close all;

% Parametrar och önskade villkor PARAMETRAR OCH ÖNSKADE VILLKOR
K1 = 0.2;                   % K1 ges, K(x) = K0 - K1*x
desired_slope = -0.51;      % Önskad sluttlutning vid x = 0.5
desired_max   = 0.255;      % Önskad maximal höjd

y0_val = 0.1;               % Randvillkor: y(0)=0.1
% Startlutningen anges initialt som s = y'(0). Vi gissar med
% exempelvis 46°: (använd tangent, eftersom lutningen = tan(vinkel))
s0_val = tand(46);          % Notera: detta ger en lutning, inte vinkelvärdet

% Parametervector: p = [K0; s]
p = [10.6753; s0_val];      % Initial gissning (från tidigare deluppgifter)

xspan = [0, 0.5];           % Intervall i x (där x mäter sidledsavstånd)
N = 1000000;                    % Antal steg i FDM (kan justeras för noggrannhet)
h = (xspan(2)-xspan(1)) / N;

tol_newton = 1e-8;          % Tolerans för Newtons metod
max_iter = 50;              % Max antal Newton-iterationer
eps_fd = 1e-6;              % Finita differenssteget för Jacobianen

% Newtons metod med finita differenser (FDM)
for iter = 1:max_iter
    % Lös IVP med FDM för nuvarande parametrar p = [K0; s]
    [x_vals, y_vals, v_vals] = fdm_solver(p, y0_val, xspan(2), N, K1);
    
    % Residual 1: Skillnaden mellan beräknad sluttlutning vid x = 0.5 och önskad lutning.
    R1 = v_vals(end) - desired_slope;
    
    % Residual 2: Skillnaden mellan den beräknade maxhöjden (över intervallet) och önskad maxhöjd.
    R2 = max(y_vals) - desired_max;
    
    F_val = [R1; R2];
    
    % Avbryt om residualerna är små nog
    if norm(F_val, inf) < tol_newton
        break;
    end
    
    % räkna jacobianen m finita differenser
    J = zeros(2,2); % tom matris för jacobian
    for i = 1:2
        dp = zeros(2,1); % Skapar en nollvektor med samma dimension som p.
        dp(i) = eps_fd;   % Sätter en liten perturbation (eps_fd) i den i:te komponenten.
        p_plus = p + dp;   % Skapar en ny parametervektor med en liten förändring i den i:te komponenten.

        % Löser ODE-systemet med den perturberade parametervektorn p_plus
        [~, y_vals_p, v_vals_p] = fdm_solver(p_plus, y0_val, xspan(2), N, K1);

        % Beräknar residualerna för den perturberade lösningen:
        % R1_p är skillnaden mellan sluttlutningen vid x=x_end och önskad sluttlutning.
        R1_p = v_vals_p(end) - desired_slope;
        
        % R2_p är skillnaden mellan den beräknade maxhöjden och den önskade maxhöjden.
        R2_p = max(y_vals_p) - desired_max;

        % Samlar de perturberade residualerna i en vektor.
        F_val_p = [R1_p; R2_p]; 
        
        % Finita differens-approximation: skillnaden (F(p+dp)-F(p)) dividerat med eps_fd 
        % ger en approximation av den partiella derivatan av F med avseende på p(i).
        % Resultatet lagras i kolumn i i Jacobianen J.
        J(:, i) = (F_val_p - F_val) / eps_fd;
    end
    
    % Newton-uppdatering: Lös J * delta = -F_val
    delta = -J \ F_val;
    p = p + delta;
end

% Extrahera de justerade parametrarna
K0_new = p(1);
s_new  = p(2);

% Lös ODE-systemet med de konvergerade parametrarna för att få hela kranens profil
[x_vals, y_vals, v_vals] = fdm_solver(p, y0_val, xspan(2), N, K1);
max_y = max(y_vals);
final_slope = v_vals(end);

% Lokal funktion: fdm_solver
function [x_vals, y_vals, v_vals] = fdm_solver(p, y0, x_end, N, K1)
    K0 = p(1);
    s = p(2);
    h = (x_end - 0) / N;
    x_vals = linspace(0, x_end, N+1)';
    
    y_vals = zeros(N+1, 1);
    v_vals = zeros(N+1, 1);
    
    % Sätt initialvillkoren
    y_vals(1) = y0;
    v_vals(1) = s;
    
    % Använd explicit Euler-metod för integration
    for i = 1:N
        x_i = x_vals(i);
        y_i = y_vals(i);
        v_i = v_vals(i);
        % Beräkna y'' enligt differentialekvationen:
        ypp = - (K0 - K1*x_i) * y_i * (1 + v_i^2)^(3/2);
        y_vals(i+1) = y_i + h * v_i;
        v_vals(i+1) = v_i + h * ypp;
    end
end

% Utskrift av resultat
disp(['s = ', num2str(s_new, '%.7f')]);
disp(['Antal Newton-iterationer: ', num2str(iter)]);
disp(['Konvergerat K0 = ', num2str(K0_new, '%.7f')]);
disp(['Konvergerad startlutning: ', num2str(rad2deg(s_new))]);

% Plottning av kranens profil (sett från sidan)
figure;
plot(x_vals, y_vals, 'b-', 'LineWidth',2);
xlabel('x [m]');
ylabel('y [m]');
title('Kranens profil (FDM)');
grid on;
