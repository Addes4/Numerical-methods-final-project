% Parametrar och startvillkor
K1           = 0.2;
Y0           = [0.1; tan(deg2rad(46))];  % [y(0); y'(0)]
xspan        = [0, 0.5];

% Målvärde för sluttlutning y'(L)
target_slope = -0.51;

% Inställningar för ODE45
opts_work = odeset('RelTol',1e-12,'AbsTol',1e-13);
opts_ref  = odeset('RelTol',1e-14,'AbsTol',1e-15);

% --- 1) Lös K0 via fzero + inline ode45/deval ---
g = @(K0) deval( ...
         ode45( ...
           @(x,y)[ ...
             y(2); ...
             -(K0 - K1*x)*y(1)*(1 + y(2)^2)^(3/2) ...
           ], ...
           xspan, Y0, opts_work ...
         ), ...
         xspan(end), 2 ...
       ) - target_slope;

K0_initial  = 11;
K0_solution = fzero(g, K0_initial);

% --- 2) Definiera ODE-systemet med funktionshandle ---
odeFun = @(x,y)[ ...
  y(2); ...
  -(K0_solution - K1*x)*y(1)*(1 + y(2)^2)^(3/2) ...
];

% --- 3) Integrera med arbets- och referensinställningar ---
[x_w, Y_w] = ode45(odeFun, xspan, Y0, opts_work);
[x_r, Y_r] = ode45(odeFun, xspan, Y0, opts_ref);

% --- 4) Rå maximihöjd från arbets- respektive referenslösning ---
[raw_max_yw, idx_w] = max(Y_w(:,1));
[raw_max_yr, idx_r] = max(Y_r(:,1));

% --- 5) 3:e ordningens interpolation kring maximalpunkt ---
if idx_w > 2 && idx_w < length(x_w)-1
    x_seg = x_w(idx_w-2:idx_w+1);
    y_seg = Y_w(idx_w-2:idx_w+1,1);
    p     = polyfit(x_seg, y_seg, 3);        % 3:e ordningens polynom
    dp    = polyder(p);                     % Derivatan är ett 2:a ordningspolynom
    r_all = roots(dp);                      % Kritiska punkter
    % Välj reella rötter inom segmentet
    r = r_all(imag(r_all)==0 & r_all>=min(x_seg) & r_all<=max(x_seg));
    if ~isempty(r)
        % Välj den som ger högst polynomvärde
        [~, imax]      = max(polyval(p, r));
        x_max_interp   = r(imax);
        max_yw_interp  = polyval(p, x_max_interp);
    else
        x_max_interp  = x_w(idx_w);
        max_yw_interp = raw_max_yw;
    end
else
    % Kantfall: använd råvärdet
    x_max_interp  = x_w(idx_w);
    max_yw_interp = raw_max_yw;
end

% --- 6) A posteriori-feluppskattning (rå) ---
method_err = abs(raw_max_yr - raw_max_yw);

% --- 7) Visa resultat ---
fprintf('Löst K0                   = %.6f\n',     K0_solution);
fprintf('Rå maximihöjd (h)          = %.6f\n',     raw_max_yw);
fprintf('Interpolerad maxhöjd (3:a) = %.6f\n',     max_yw_interp);
fprintf('Maxhöjd (h/2) (rå)         = %.6f\n',     raw_max_yr);
fprintf('Metodfel (rå)              = %.2e\n',     method_err);
