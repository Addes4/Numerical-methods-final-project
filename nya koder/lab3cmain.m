function solve_deluppgift3(K0)
    % Parametrar
    K0 = 10.675320898642
    K1 = 0.2;
    y0 = 0.1;
    s0 = tan(deg2rad(46));  % y'(0) = tan(46°)
    Y0 = [y0; s0];
    x0 = 0;
    L = 0.5;
    h = 0.00001;  % Steglängd
    
    % Anonym funktion med inbäddade parametrar
    f = @(x, Y) ode_system(x, Y, K0, K1);
    
    % Lös systemet med RK4
    [x_vals, Y] = rk4_system(f, x0, Y0, h, L);
    
    % Extrahera y-lösningarna
    y_vals = Y(1, :);
    
    % Hitta index för det största y-vä  rdet
    [~, idx] = max(y_vals);
    
    x1 = x_vals(idx-2); y1 = y_vals(idx-2);
    x2 = x_vals(idx-1); y2 = y_vals(idx-1);
    x3 = x_vals(idx);   y3 = y_vals(idx);
    x4 = x_vals(idx+1); y4 = y_vals(idx+1);
    
    % Andragradsanpassning (befintlig)
    A_quad = [x2^2, x2, 1;
              x3^2, x3, 1;
              x4^2, x4, 1];
    b_quad = [y2; y3; y4];
    coeff_quad = A_quad\b_quad;
    a_q = coeff_quad(1); b_q = coeff_quad(2);
    x_max_quad = -b_q/(2*a_q);
    
    % Tredjegradsanpassning
    A_cubic = [x1^3, x1^2, x1, 1;
               x2^3, x2^2, x2, 1;
               x3^3, x3^2, x3, 1;
               x4^3, x4^2, x4, 1];
    b_cubic = [y1; y2; y3; y4];
    coeff_cubic = A_cubic\b_cubic;
    a_c = coeff_cubic(1); b_c = coeff_cubic(2); c_c = coeff_cubic(3);
    % Derivata av kubiskt polynom: 3a_c*x² + 2b_c*x + c_c
    x_max_cubic = roots([3*a_c, 2*b_c, c_c]);
    x_max_cubic = x_max_cubic(imag(x_max_cubic)==0); % Välj riktigt rot

% Feluppskattning
error_x = abs(x_max_quad - x_max_cubic);
fprintf('Feluppskattning: %.3e\n', error_x);
