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
    
    % Hitta index för det största y-värdet
    [~, idx] = max(y_vals);
    
    % Välj tre punkter kring maximum: en före, vid och efter maximum
    x1 = x_vals(idx-1); y1 = y_vals(idx-1);
    x2 = x_vals(idx);   y2 = y_vals(idx);
    x3 = x_vals(idx+1); y3 = y_vals(idx+1);
    
    % Sätt upp systemet för ett andragradspolynom: y = a*x^2 + b*x + c
    A = [x1^2, x1, 1;
         x2^2, x2, 1;
         x3^2, x3, 1];
    b_vec = [y1; y2; y3];
    
    % Lös för koefficienterna a, b och c
    coeff = A\b_vec;
    a = coeff(1); 
    b = coeff(2); 
    c = coeff(3);
    
    % Beräkna vertex för polynomet: x_max = -b/(2a)
    x_max = -b / (2*a);
    y_max = a*x_max^2 + b*x_max + c;
    
    % Skriv ut resultatet
    fprintf('Interpolerat maximum (deluppgift c):\n');
    fprintf('  x = %.9f\n', x_max);
    fprintf('  y = %.6f\n', y_max);
end
