function lab3a_adv
    % Parametrar
    K0 = 11;
    K1 = 0.2;

    % Initialvillkor
    y0 = 0.1;
    s0 = tan(deg2rad(46));  % y'(0)
    
    % Intervall
    xspan = [0 0.5];

    % Avancerad: använd ode45 med två toleransnivåer för metodfel
    opts1 = odeset('RelTol',1e-13,'AbsTol',1e-14);
    [~,Y1] = ode45(@(x,Y) ode_system(x,Y,K0,K1), xspan, [y0; s0], opts1);
    y1 = Y1(end,1);
    
    opts2 = odeset('RelTol',1e-15,'AbsTol',1e-16);
    [~,Y2] = ode45(@(x,Y) ode_system(x,Y,K0,K1), xspan, [y0; s0], opts2);
    y2 = Y2(end,1);
    
    method_err = abs(y2 - y1);

    % Indatafel ±1%
    params_nom    = [y0, s0, K0, K1];
    err_inputs    = zeros(1,4);
    for j = 1:4
        p = params_nom;
        p(j) = p(j) * 1.01;
        [~,Yp] = ode45(@(x,Y) ode_system(x,Y,p(3),p(4)), xspan, [p(1); p(2)], opts1);
        err_inputs(j) = abs(Yp(end,1) - y2);
    end
    input_err = sum(err_inputs);
    
    total_err = method_err + input_err;

    % Utskrift
disp('Avancerad deluppgift a)')
disp(['  y(0.5)       = ' num2str(y1,        '%.8f') ' m'])
disp(['  Metodfel     = ' num2str(method_err, '%.2e')  ' m'])
disp(['  Indatafel ≤=  ' num2str(input_err,  '%.2e')  ' m'])
disp(['  Total fel ≤=  ' num2str(total_err,  '%.2e')  ' m'])
end


function dYdx = ode_system(x, Y, K0, K1)
    y = Y(1);
    v = Y(2);
    dYdx = [v; - (K0 - K1*x) * y * (1 + v^2)^(3/2)];
end
