clearvars
sim.GPUAccel = 1;

% Physical constants

ech = 1.602176634e-19;  % electron charge, C
amu = 1.66053906660e-27;    % atomic mass unit, kg
eps0 = 8.8541878128e-12;    % vacuum electric permittivity

% Trap parameters

Vac = 700;  % AC-voltage, V
Udc = 2;   % DC radial voltage, V
RF = 25e6; % Radiofrequency, Hz

r0 = 0.624e-3; % radial geometric constant, m
z0 = 2.25e-3; % axial geometric constant, m
kap = 0.0567;    % geometric factor
asym = 1;   % radial asymmetry factor

ax = -4*ech*kap*Udc/(40*amu)/z0^2/(RF*2*pi)^2;
ay = ax;
az = -2*ax;
qx = 2*ech*Vac/(40*amu)/r0^2/(RF*2*pi)^2;
qy = -qx;

x_eq=[];
y_eq=[];
z_eq=[];
norm_modes = [];
frs = [];
w_n = [];
l = [];
Ca_40_ind = [];
AncillaIon_ind = [];

masses = [40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40];
chars = ones(1, size(masses, 2));
N = size(masses, 2);

[x_eq, y_eq, z_eq, norm_modes, frs, w_n, l] = get_modes(masses, chars, RF, ax, qx);

[xData, yData] = prepareCurveData( [], z_eq );
ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.234419756878536 0.206164801457355];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
fitresult.a
%%
average_distance = [1.8593e-05, ];