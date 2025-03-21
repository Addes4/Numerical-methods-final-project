function faucet_profile()
    % Parametrar
    K0 = 11;
    K1 = 0.2;

    % Startvillkor
    y0 = 0.1;                   % Handfatets kant (i meter)
    theta0 = 46 * pi/180;       % Startvinkel i radianer (46°)
    yp0 = tan(theta0);          % y'(0) = tan(theta0)

    % ODE-system: u(1) = y, u(2) = y'
    odefun = @(x, u) [ u(2);
                      -(K0 - K1*x) * u(1) * (1 + u(2)^2)^(3/2) ];

    % Lös ODE:n med ode45
    [x, u] = ode45(odefun, [0, 0.5], [y0; yp0]);

    % Sluthöjd vid x = 0.5
    y_end = u(end, 1);
    fprintf('Vid x = 0.5 avslutas kranen på höjden y = %.8f m\n', y_end);

    % Plotta profilen
    figure;
    plot(x, u(:,1), 'b-', 'LineWidth', 2);
    xlabel('x (m)');
    ylabel('y (m)');
    title('Profil av den böjda vattenkranen');
    grid on;
end