%% Initializing constants
clearvars
sim.GPUAccel = 1;

% Physical constants

ech = 1.602176634e-19;  % electron charge, C
amu = 1.66053906660e-27;    % atomic mass unit, kg
eps0 = 8.8541878128e-12;    % vacuum electric permittivity
pi = 3.14159265359;

% Trap parameters

Vac = 850;  % AC-voltage, V
Udc = 3;   % DC radial voltage, V
RF = 25e6; % Radiofrequency, Hz

r0 = 0.624e-3; % radial geometric constant, m
z0 = 2.25e-4; % axial geometric constant, m
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

%% Initializung chain parameters
masses = [40 40 48 48 48 40 40 40 48 48];
tweez_w = [0 0 0 0 0 ];
chars = ones(1, size(masses, 2));
    
N = size(masses, 2);    

[x, y, z, x_eq, y_eq, z_eq, norm_modes, frs, w_n, l, time] = get_modes(masses, chars ,RF, ax, qx); 
heatmap_norm_modes(masses, norm_modes, frs, w_n, x_eq, y_eq, z_eq);
%%
axial_incr_ancilla48 = [1.1939 1.3184 1.4218 1.5148 1.6013 1.6833 1.7621 1.8383 1.9124 1.9848];
axial_incr_in = [2.0223 2.0943 2.1461 2.1919 2.2343 2.2744];
periodic = [2.4540 2.5926 2.7039 2.8021 2.8910 2.9728];
plot(axial_incr_ancilla48*w_n, 'bo-');
title('Axial COM mode increasing');
xlabel('Number of ions in each ancilla group');
ylabel('Axial COM mode, Hz');
grid on;
hold on;
plot(axial_incr_in*w_n, 'ro-')
plot(periodic*w_n, 'ko-')
hold off;
legend('edge ancillas', 'inner ancillas (3 groups)', 'inner ancillas but with smaller period (7 groups)')