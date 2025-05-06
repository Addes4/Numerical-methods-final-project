% parametrar
K0 = 11;
K1 = 0.2;
y0 = 0.1;
s0 = tan(deg2rad(46));
Y0 = [y0; s0];

% steglängd o int
x0 = 0;
L = 0.5;
h = 1e-3;
h2 = h/2;
h4 = h / 4;


%  anropa ode
f = @(x,Y) ode_system(x, Y, K0, K1);

% lös m olika steglängd 
[~, Yh] = rk4_system(f, x0, Y0, h,  L);
[~, Yh2] = rk4_system(f, x0, Y0, h2, L);
[~, Yh4] = rk4_system(f, x0, Y0, h4, L);

yh = Yh(1,end);   % y vid x = 0.5 för steglängd h
yh2 = Yh2(1,end); 
yh4  = Yh4(1,end); 

% fel mellan lösningarna 
E1 = abs(yh  - yh2);
E2 = abs(yh2 - yh4);

% konv.ordn
p_est = log(E1/E2) / log(2);

disp(['y(0.5)  = ' num2str(yh2, '%.8f')])
disp(['metodfelet = ' num2str(E1,  '%.2e')])
disp([' p  = ' num2str(p_est, '%.4f')])
