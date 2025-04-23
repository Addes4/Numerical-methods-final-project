function deluppgift_d_builtin()
% Basnivå med inbyggda ode45 & fzero, inkl. metodfelsskattning

%% Parametrar och mål
x0            = 0;                      % Startpunkt x
L             = 0.5;                    % Slutpunkt x
y0            = 0.1;                    % Starthöjd y(0)
K1            = 0.2;                    % Konstant i krökningsfunktion

target_slope  = -0.51;                  % Målvärde y'(L)
target_height = 0.255;                  % Målvärde maxhöjd y

% Toleranser för metodfelsskattning
optsWork = odeset('RelTol',1e-6,'AbsTol',1e-8);
optsRef  = odeset('RelTol',1e-9,'AbsTol',1e-10);

% 1) Arbetslösning
% a) Hitta s0 så att max y = target_height
t_funS0 = @(s0) getMaxY(findK0(s0, optsWork, x0, L, y0, K1), s0, optsWork, x0, L, y0, K1) - target_height;
s0_work  = fzero(t_funS0, tan(deg2rad(46)));

% b) Hitta K0 så att y'(L) = target_slope
k_funK0 = @(K0) getSlope(K0, s0_work, optsWork, x0, L, y0, K1) - target_slope;
K0_work  = fzero(k_funK0, 11);

% c) Simulera slutvärden
dataWork    = simulate(K0_work, s0_work, optsWork, x0, L, y0, K1);
slope_work  = dataWork.slope;
max_work    = dataWork.maxY;

% 2) Referenslösning
r_funS0 = @(s0) getMaxY(findK0(s0, optsRef, x0, L, y0, K1), s0, optsRef, x0, L, y0, K1) - target_height;
s0_ref  = fzero(r_funS0, s0_work);
r_funK0 = @(K0) getSlope(K0, s0_ref, optsRef, x0, L, y0, K1) - target_slope;
K0_ref  = fzero(r_funK0, K0_work);
dataRef    = simulate(K0_ref, s0_ref, optsRef, x0, L, y0, K1);
slope_ref  = dataRef.slope;
max_ref    = dataRef.maxY;

% 3) Metodfelsskattning
method_err_slope  = abs(slope_ref - slope_work);
method_err_height = abs(max_ref   - max_work);

% 4) Utskrift
disp('Deluppgift d) med ode45 & fzero')
disp(['  Startlutning s0    = ' num2str(s0_work,      '%.6f')])
disp(['  Parameter K0       = ' num2str(K0_work,      '%.6f')])
disp(['  y''(L) (arb)       = ' num2str(slope_work,   '%.6f')])
disp(['  max y (arb)        = ' num2str(max_work,     '%.6f') ' m'])
disp(['  Metodfel y''       = ' num2str(method_err_slope,  '%.2e')])
disp(['  Metodfel max y     = ' num2str(method_err_height,'%.2e')])
disp(['  Startlutning s0    = ' num2str(rad2deg(s0_work),      '%.6f')])

% Lokala funktioner
    function out = simulate(K0, s0, opts, x0, L, y0, K1)
        odeFun = @(x, Y) [Y(2); -(K0 - K1*x)*Y(1)*(1 + Y(2)^2)^(3/2)];
        [~, Y]  = ode45(odeFun, [x0 L], [y0; s0], opts);
        out.slope = Y(end,2);
        out.maxY  = max(Y(:,1));
    end

    function v = getSlope(K0, s0, opts, x0, L, y0, K1)
        v = simulate(K0, s0, opts, x0, L, y0, K1).slope;
    end

    function v = getMaxY(K0, s0, opts, x0, L, y0, K1)
        v = simulate(K0, s0, opts, x0, L, y0, K1).maxY;
    end

    function K0out = findK0(s0, opts, x0, L, y0, K1)
        funK = @(K0) getSlope(K0, s0, opts, x0, L, y0, K1) - target_slope;
        K0out = fzero(funK, 11);
    end
end
