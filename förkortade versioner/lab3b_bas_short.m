
% parametrra
K1 = 0.2;
x0 = 0;
L = 0.5;
Y0 = [0.1; tan(deg2rad(46))];
target_slope = -0.51;
h = 1e-3; % steg

% sekant gissningar
K0_prev = 9;
K0_curr = 11;

% 1a funktionsvärden (slutlutn - target)
f_prev = compute_metrics(K0_prev, Y0(2), x0, Y0, h, L, K1) - target_slope;
f_curr = compute_metrics(K0_curr, Y0(2), x0, Y0, h, L, K1) - target_slope;

% tol o iter
tol = 1e-5;
iter = 0;
maxIter = 50;

% konvergensanalys
error_list = [];

% sekantloop
while abs(K0_curr - K0_prev) > tol && iter < maxIter
    K0_new = K0_curr - f_curr * (K0_curr - K0_prev) / (f_curr - f_prev);

    % spara fel för konv.ord
    error_list(end+1) = abs(K0_curr - K0_prev);

    % upd
    K0_prev = K0_curr;
    f_prev = f_curr;
    K0_curr = K0_new;
    f_curr = compute_metrics(K0_curr, Y0(2), x0, Y0, h, L, K1) - target_slope;

    iter = iter + 1;
end

K0_solution = K0_curr;

%metodfel
y_work = compute_metrics(K0_solution, Y0(2), x0, Y0, h,   L, K1);
y_ref  = compute_metrics(K0_solution, Y0(2), x0, Y0, h/2, L, K1);
method_err = abs(y_ref - y_work);

% konv.ord beräkning
if length(error_list) >= 3
    e_n   = error_list(end);
    e_nm1 = error_list(end-1);
    e_nm2 = error_list(end-2);
    p = log(e_n / e_nm1) / log(e_nm1 / e_nm2);
else
    p = NaN;
end

disp(['K0 = ', num2str(K0_solution, '%.6f')])
disp(['y''(L) med h  = ', num2str(y_work, '%.6f')])
disp(['y''(L) med h/2 = ', num2str(y_ref, '%.6f')])
disp(['metodfelet = ', num2str(method_err, '%.2e')])
disp(['antal iterationer: ', num2str(iter)])
disp(['konvergensordning (p) ≈ ', num2str(p, '%.4f')])
