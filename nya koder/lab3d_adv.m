% Parametrar
x0           = 0;
L            = 0.5;
y0           = 0.1;
K1           = 0.2;
target_slope = -0.51;
target_height= 0.255;

% ----- Enradiga hjälpfunktioner -----
% Hjälpfunktion för att beräkna sluttlutning y'(L)
get_slope = @(K0,s0) ... % Anonym funktionshandle
    deval( ...   % Evaluera lösningen vid x = L
      ode45( ...  % Använd MATLAB:s ode45
      @(x,Y)[Y(2);... % Definiera ODE: Y(1)=y, Y(2)=y'
      -(K0-K1*x)*Y(1)*(1+Y(2)^2)^(3/2)], ... 
             [x0, L], [y0; s0] ), ... % Tidsintervall och startvärden
      L, 2 ); % Välj komponent 2 = y'(L)

% Hjälpfunktion för att beräkna maxhöjd y_{max}
get_height = @(K0,s0) ... % Anonym funktionshandle
    max( ...  % Hitta största y-värde
      deval( ...  % Evaluera lösningen på ett finmaskigt nät
        ode45( @(x,Y)[Y(2); -(K0-K1*x)*Y(1)*(1+Y(2)^2)^(3/2)], ...
               [x0, L], [y0; s0] ), ...
        linspace(x0, L, 500), ... % 500 punkter mellan x0 och L
        1 ) ); %Välj komponent 1 = y(x)

% ----- Rotsökning -----
% 1) Hitta startlutning s0 så att maxhöjden blir rätt
find_K0 = @(s0)... % K0 som funktion av s0
    fzero(...  % Hitta rot för sluttlutningsfelet
    @(K0) get_slope(K0,s0) - target_slope, 10);
height_err = @(s0)... % Höjdfel som funktion av s0
    get_height(find_K0(s0), s0) - target_height;
s0_sol = fzero(... % Rot för höjdfunktionen
    height_err, [tan(deg2rad(40)), tan(deg2rad(50))]);

% 2) Hitta K0 med funnet s0
K0_sol = find_K0(s0_sol); % K0 som ger rätt sluttlutning

fprintf('Resultat: s0 = %.4f°  ,  K0 = %.6f\n', rad2deg(s0_sol), K0_sol);

% ----- Beräkning av metodfel -----
% Två uppsättningar toleranser för ode45
opts1 = odeset('RelTol',1e-11,  'AbsTol',1e-12);
opts2 = odeset('RelTol',1e-13, 'AbsTol',1e-14);

% Lagra ODE-funktionen
odeFun = @(x,Y,K0) [Y(2); -(K0-K1*x)*Y(1)*(1+Y(2)^2)^(3/2)];

% 1) Med första toleransen
sol1 = ode45(@(x,Y) odeFun(x,Y,K0_sol), [x0, L], [y0; tan(s0_sol)], opts1);
s1   = deval(sol1, L, 2);
Yv1  = deval(sol1, linspace(x0, L, 500));
h1   = max(Yv1(1,:));

% 2) Med andra (striktare) toleransen
sol2 = ode45(@(x,Y) odeFun(x,Y,K0_sol), [x0, L], [y0; tan(s0_sol)], opts2);
s2   = deval(sol2, L, 2);
Yv2  = deval(sol2, linspace(x0, L, 500));
h2   = max(Yv2(1,:));

% Metodfel som skillnad mellan resultaten
err_slope  = abs(s2 - s1);
err_height = abs(h2 - h1);

fprintf('Metodfel: err_slope = %.2e, err_height = %.2e\n', err_slope, err_height);
