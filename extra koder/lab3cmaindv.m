function solve_deluppgift3(K0)
    % Parametrar
    K1 = 0.2;
    y0 = 0.1;
    s0 = tan(deg2rad(46));  % y'(0) = tan(46°)
    Y0 = [y0; s0];
    x0 = 0;
    L = 0.5;
    h = 0.00001;  % Steglängd

    % ODE-systemfunktion
    f = @(x, Y) ode_system(x, Y, K0, K1);

    % Lös systemet med RK4
    [x_vals, Y] = rk4_system(f, x0, Y0, h, L);

    % Extrahera y-lösningarna
    y_vals = Y(1, :);

    % Hitta index för största y-värdet
    [~, idx] = max(y_vals);

    % Hämta fyra punkter kring maximum
    x1 = x_vals(idx-2); y1 = y_vals(idx-2);
    x2 = x_vals(idx-1); y2 = y_vals(idx-1);
    x3 = x_vals(idx);   y3 = y_vals(idx);
    x4 = x_vals(idx+1); y4 = y_vals(idx+1);

    % Andragradsanpassning med polyfit
    x_quad = [x2, x3, x4];
    y_quad = [y2, y3, y4];
    coeff_quad = polyfit(x_quad, y_quad, 2);
    a_q = coeff_quad(1); b_q = coeff_quad(2);
    x_max_quad = -b_q / (2 * a_q);

    % Tredjegradsanpassning med polyfit
    x_cubic = [x1, x2, x3, x4];
    y_cubic = [y1, y2, y3, y4];
    coeff_cubic = polyfit(x_cubic, y_cubic, 3);

    % Derivera tredjegradspolynom
    deriv_coeffs = polyder(coeff_cubic);
    x_roots = roots(deriv_coeffs);
    x_roots = x_roots(imag(x_roots) == 0);  % Endast reella rötter

    % Välj maximum av polynomet bland reella rötter
    y_root_vals = polyval(coeff_cubic, x_roots);
    [~, i_max] = max(y_root_vals);
    x_max_cubic = x_roots(i_max);

    % Feluppskattning
    error_x = abs(x_max_quad - x_max_cubic);
    fprintf('Feluppskattning: %.3e\n', error_x);
end

