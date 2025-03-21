% Water Tap Problem

% Part a: Solve IVP and find y(0.5)
K0_a = 11;
K1 = 0.2;
s0_a = tan(46 * pi/180);
xspan = [0, 0.5];
[y_end_a, yp_end_a, ~] = solve_ode(K0_a, s0_a, K1, xspan);
fprintf('Part a: y(0.5) = %.6f\n', y_end_a);

% Part b: Adjust K0 to get y''(0.5) = -0.51
target_yp_end_b = -0.51;
epsilon_b = 0.01;
tolerance_b = 1e-6;
max_iter_b = 100;
K0_initial_b = K0_a;

K0_current = K0_initial_b;
for iter = 1:max_iter_b
    [~, yp_end, ~] = solve_ode(K0_current, s0_a, K1, xspan);
    F_current = yp_end - target_yp_end_b;
    if abs(F_current) < tolerance_b
        break;
    end
    % Compute F at K0 + epsilon
    [~, yp_plus, ~] = solve_ode(K0_current + epsilon_b, s0_a, K1, xspan);
    F_plus = yp_plus - target_yp_end_b;
    dFdK0 = (F_plus - F_current) / epsilon_b;
    % Update K0
    delta_K0 = -F_current / dFdK0;
    K0_current = K0_current + delta_K0;
end
K0_b = K0_current;
fprintf('Part b: K0 = %.6f\n', K0_b);

% Part c: Find maximum height with new K0
[~, ~, y_max_c] = solve_ode(K0_b, s0_a, K1, xspan);
fprintf('Part c: Maximum height y_max = %.6f\n', y_max_c);

% Part d: Adjust K0 and s0 to meet y_max = 0.255 and y''(0.5) = -0.51
target_yp_end_d = -0.51;
target_y_max_d = 0.255;
epsilon_d = 0.01;
tolerance_d = 1e-6;
max_iter_d = 20;

% Initial guess: K0 from part b and s0 from part a
X = [K0_b; s0_a];

for iter = 1:max_iter_d
    [~, yp_end, y_max] = solve_ode(X(1), X(2), K1, xspan);
    F1 = yp_end - target_yp_end_d;
    F2 = y_max - target_y_max_d;
    F = [F1; F2];
    if norm(F) < tolerance_d
        break;
    end
    % Compute Jacobian
    % Perturb K0
    [~, yp_K_plus, y_max_K_plus] = solve_ode(X(1) + epsilon_d, X(2), K1, xspan);
    F1_K_plus = yp_K_plus - target_yp_end_d;
    F2_K_plus = y_max_K_plus - target_y_max_d;
    dF1_dK0 = (F1_K_plus - F1) / epsilon_d;
    dF2_dK0 = (F2_K_plus - F2) / epsilon_d;
    
    % Perturb s0
    [~, yp_s_plus, y_max_s_plus] = solve_ode(X(1), X(2) + epsilon_d, K1, xspan);
    F1_s_plus = yp_s_plus - target_yp_end_d;
    F2_s_plus = y_max_s_plus - target_y_max_d;
    dF1_ds0 = (F1_s_plus - F1) / epsilon_d;
    dF2_ds0 = (F2_s_plus - F2) / epsilon_d;
    
    J = [dF1_dK0, dF1_ds0;
         dF2_dK0, dF2_ds0];
    
    % Solve for delta
    delta = J \ (-F);
    
    % Update parameters
    X = X + delta;
end

K0_d = X(1);
s0_d = X(2);
fprintf('Part d: K0 = %.6f, s0 = %.6f (initial angle = %.2f degrees)\n', ...
        K0_d, s0_d, atan(s0_d)*180/pi);

% Helper function to solve the ODE
function [y_end, yp_end, y_max] = solve_ode(K0, s0, K1, xspan)
    y0 = [0.1; s0];
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
    sol = ode45(@(x, u) odefun(x, u, K0, K1), xspan, y0, options);
    % Ensure solution reached x=0.5
    if sol.x(end) < xspan(2)
        error('ODE solver did not reach x=0.5');
    end
    u_end = deval(sol, xspan(2));
    y_end = u_end(1);
    yp_end = u_end(2);
    % Find maximum y
    x_sol = sol.x;
    y_sol = sol.y(1,:);
    yp_sol = sol.y(2,:);
    % Detect sign change in yp
    sign_changes = find(diff(sign(yp_sol)) < 0);
    if ~isempty(sign_changes)
        idx = sign_changes(1);
        x1 = x_sol(idx);
        x2 = x_sol(idx+1);
        yp1 = yp_sol(idx);
        yp2 = yp_sol(idx+1);
        % Linear interpolation for x where yp=0
        x_max = x1 - yp1 * (x2 - x1) / (yp2 - yp1);
        y_max = deval(sol, x_max);
        y_max = y_max(1);
    else
        y_max = max(y_sol);
    end
end

% ODE definition
function dudx = odefun(x, u, K0, K1)
    y = u(1);
    yp = u(2);
    K = K0 - K1 * x;
    ypp = -K * y * (1 + yp^2)^(3/2);
    dudx = [yp; ypp];
end