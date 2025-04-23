% lab3alltindata.m: RK4-kran med felskattning och störningsanalys
% =================================================================================
% Script som löser ODE-system för vattenkran, skattar metodfel, och gör störningsanalys.

%% Huvudprogram
% Basdata
data_base.K0            = 11;
data_base.K1            = 0.2;
data_base.y0            = 0.1;
data_base.s0            = tan(deg2rad(46));
data_base.target_slope  = -0.51;
data_base.target_height = 0.255;

% Kör utan störning
fprintf('--- Resultat utan störning ---\n');
res_base = run_all(data_base);

% Lägg på 1%% störning på all indata
data_pert = data_base;
fields = fieldnames(data_base);
for i = 1:numel(fields)
    data_pert.(fields{i}) = data_base.(fields{i}) * 1.01;
end

fprintf('--- Resultat med 1%% störning ---\n');
res_pert = run_all(data_pert);

% Jämför och visa differenser
fprintf('--- Störningsdifferenser ---\n');
% absoluta differenser
d_y   = abs(res_pert.a_y    - res_base.a_y);
d_K0b = abs(res_pert.b_K0   - res_base.b_K0);
d_K0c = abs(res_pert.c_K0   - res_base.c_K0);
d_s0  = abs(res_pert.d_s0K0(1) - res_base.d_s0K0(1));
d_K0d = abs(res_pert.d_s0K0(2) - res_base.d_s0K0(2));
% procentuella förändringar
p_y   = d_y   / abs(res_base.a_y)      * 100;
p_K0b = d_K0b / abs(res_base.b_K0)     * 100;
p_K0c = d_K0c / abs(res_base.c_K0)     * 100;
p_s0  = d_s0  / abs(res_base.d_s0K0(1)) * 100;
p_K0d = d_K0d / abs(res_base.d_s0K0(2)) * 100;

fprintf('a) Δy(0.5)      = %.2e (%.2f%%)\n', d_y,   p_y);
fprintf('b) ΔK0          = %.2e (%.2f%%)\n', d_K0b, p_K0b);
fprintf('c) ΔK0 (del c)  = %.2e (%.2f%%)\n', d_K0c, p_K0c);
fprintf('d) Δs0          = %.2e (%.2f%%)\n', d_s0,  p_s0);
fprintf('d) ΔK0 (del d)  = %.2e (%.2f%%)\n', d_K0d, p_K0d);

%% Lokala funktioner
function res = run_all(data)
    res.a_y    = run_sub_a(data);
    res.b_K0   = run_sub_b(data);
    res.c_K0   = run_sub_c(data);
    res.d_s0K0 = run_sub_d(data);
end

function y_end = run_sub_a(data)
    h  = 1e-4; h2 = h/2;
    x0 = 0; L = 0.5;
    Y0 = [data.y0; data.s0];
    f = @(x,Y) crane_ode(x,Y,data.K0,data.K1);
    [~,Yh ] = rk4_system(f,x0,Y0,h, L);
    [~,Yh2] = rk4_system(f,x0,Y0,h2,L);
    y_end = Yh2(1,end);
end

function K0_sol = run_sub_b(data)
    x0 = 0; L = 0.5; h = 1e-3;
    Y0 = [data.y0; data.s0];
    K0_prev = data.K0 * 0.5;
    K0_curr = data.K0 * 1.5;
    f_prev = compute_end_slope(K0_prev, x0, Y0, h, L, data.K1) - data.target_slope;
    f_curr = compute_end_slope(K0_curr, x0, Y0, h, L, data.K1) - data.target_slope;
    tol = 1e-5; iter = 0; maxI = 50;
    while abs(K0_curr - K0_prev) > tol && iter < maxI
        K0_new  = K0_curr - f_curr*(K0_curr - K0_prev)/(f_curr - f_prev);
        K0_prev = K0_curr; f_prev = f_curr; K0_curr = K0_new;
        f_curr  = compute_end_slope(K0_curr, x0, Y0, h, L, data.K1) - data.target_slope;
        iter = iter + 1;
    end
    K0_sol = K0_curr;
end

function K0_sol = run_sub_c(data)
    x0 = 0; L = 0.5; h = 1e-3;
    Y0 = [data.y0; data.s0];
    K0_prev = data.K0 * 0.5;
    K0_curr = data.K0 * 1.5;
    f_prev  = compute_slope(K0_prev, x0, Y0, h, L, data.K1) - data.target_slope;
    f_curr  = compute_slope(K0_curr, x0, Y0, h, L, data.K1) - data.target_slope;
    tol = 1e-4; iter = 0; maxI = 50;
    while abs(K0_curr - K0_prev) > tol && iter < maxI
        K0_new  = K0_curr - f_curr*(K0_curr - K0_prev)/(f_curr - f_prev);
        K0_prev = K0_curr; f_prev = f_curr; K0_curr = K0_new;
        f_curr  = compute_slope(K0_curr, x0, Y0, h, L, data.K1) - data.target_slope;
        iter = iter + 1;
    end
    K0_sol = K0_curr;
end

function out = run_sub_d(data)
    x0 = 0; L = 0.5; h = 1e-5;
    s0_prev = data.s0 * 0.8;
    s0_curr = data.s0 * 1.2;
    tol_s = 1e-5; iter = 0; maxI = 50;
    while abs(s0_curr - s0_prev) > tol_s && iter < maxI
        K0_p = find_K0_for_s0(s0_prev, data.target_slope, x0, [data.y0; 0], h, L, data.K1);
        [~,h_p] = compute_metrics(K0_p, s0_prev, x0, [data.y0; 0], h, L, data.K1);
        K0_c = find_K0_for_s0(s0_curr, data.target_slope, x0, [data.y0; 0], h, L, data.K1);
        [~,h_c] = compute_metrics(K0_c, s0_curr, x0, [data.y0; 0], h, L, data.K1);
        f_p = h_p - data.target_height; f_c = h_c - data.target_height;
        s_new = s0_curr - f_c*(s0_curr - s0_prev)/(f_c - f_p);
        s0_prev = s0_curr; s0_curr = s_new; iter = iter + 1;
    end
    s0_sol = s0_curr;
    K0_sol = find_K0_for_s0(s0_sol, data.target_slope, x0, [data.y0; 0], h, L, data.K1);
    out = [s0_sol, K0_sol];
end

function dYdx = crane_ode(x, Y, K0, K1)
    y = Y(1); v = Y(2);
    dYdx = [v; -(K0 - K1*x)*y*(1+v^2)^(3/2)];
end

function [x, Y] = rk4_system(f, x0, Y0, h, L)
    x = x0:h:L; N = numel(x);
    Y = zeros(length(Y0), N); Y(:,1) = Y0;
    for i=1:N-1
        k1=f(x(i),Y(:,i)); k2=f(x(i)+h/2, Y(:,i)+h*k1/2);
        k3=f(x(i)+h/2, Y(:,i)+h*k2/2); k4=f(x(i)+h, Y(:,i)+h*k3);
        Y(:,i+1)=Y(:,i)+(h/6)*(k1+2*k2+2*k3+k4);
    end
end

function s_end = compute_end_slope(K0, x0, Y0, h, L, K1)
    f=@(x,Y) crane_ode(x,Y,K0,K1);
    [~,Y]=rk4_system(f,x0,Y0,h,L); s_end=Y(2,end);
end

function slope_end = compute_slope(K0, x0, Y0, h, L, K1)
    slope_end=compute_end_slope(K0,x0,Y0,h,L,K1);
end

function [s_end, h_max] = compute_metrics(K0, s0, x0, Y0, h, L, K1)
    Y0_loc=[Y0(1);s0]; f=@(x,Y) crane_ode(x,Y,K0,K1);
    [~,Y]=rk4_system(f,x0,Y0_loc,h,L); s_end=Y(2,end); h_max=max(Y(1,:));
end

function K0_out = find_K0_for_s0(s0, target_slope, x0, Y0, h, L, K1)
    K0_prev=5; K0_curr=15; tol=1e-4; iter=0; maxI=50;
    [s_p,~]=compute_metrics(K0_prev,s0,x0,Y0,h,L,K1); f_p=s_p-target_slope;
    [s_c,~]=compute_metrics(K0_curr,s0,x0,Y0,h,L,K1); f_c=s_c-target_slope;
    while abs(K0_curr-K0_prev)>tol && iter<maxI
        K0_new=K0_curr-f_c*(K0_curr-K0_prev)/(f_c-f_p);
        K0_prev=K0_curr; f_p=f_c; K0_curr=K0_new;
        [s_c,~]=compute_metrics(K0_curr,s0,x0,Y0,h,L,K1); f_c=s_c-target_slope;
        iter=iter+1;
    end
    K0_out=K0_curr;
end
