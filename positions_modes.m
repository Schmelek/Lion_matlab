
%% General code calculating the mixed species coulomb crystal normal modes
% This script creates a multi species ion coulomb crystal and calculates 
% equilibrium positions coordinates of ions. Those equilibrium positions
% are used for calculation of Hessian matrix, which eigenvalues give the
% normal mode spectrum. The modes are calculated with respect to given ion's
% secular frequency.

% clear the workspace

clearvars

% Physical constants

ech = 1.602176634e-19;  % electron charge, C
amu = 1.66053906660e-27;    % atomic mass unit, kg
eps0 = 8.8541878128e-12;    % vacuum electric permittivity

%% If you want to specify trap parameters fully, do it below. 
% Parameters a and q are calculated for Ca(40) ion.

% RF = 5e6; % Radiofrequency, Hz
% r0 = 1e-3; % radial geometric constant, m
% z0 = 1e-3; % axial geometric constant, m
% kap = 0.244;    % geometric factor
% asym = 1;   % radial asymmetry factor
% Vac = 120;  % AC-voltage, V/m
% Udc = 60;   % DC radial voltage, V/m
% Uz = 20;   % DC axial voltage, V/m
% 
% ax = 4*ech*Udc/(40*amu)/r0^2/(RF*2*pi)^2-4*ech*kap*Uz/(40*amu)/z0^2/(RF*2*pi)^2;
% ay = -4*ech*Udc/(40*amu)/r0^2/(RF*2*pi)^2-4*ech*kap*Uz/(40*amu)/z0^2/(RF*2*pi)^2;
% az = 8*ech*kap*Uz/(40*amu)/z0^2/(RF*2*pi)^2;
% qx = 2*ech*Vac/(40*amu)/r0^2/(RF*2*pi)^2;
% qy = -qx;

%% If you want to specify only symmetric (a,q) Mathieu parameters, do it here.
% choose it if the actual trap parameters are not needed

RF = 5.634e6; % Radiofrequency, Hz
qx = 0.03;  % q Mathieu parameter for Ca(40) or the first ion
ax = -0.00001;  % a Mathieu parameter for Ca(40) or the first ion
qy = qx;
ay = ax;
az = -2*ax;
%% If you want to specify secular frequencies directly (in Hz!!!), do it here.

% RF = 5e6; % Radiofrequency, Hz
% wx = 12.26e6; % x-axis secular frequency, Hz
% wy = 11.19e6; % y-axis secular frequency, Hz
% wz = 2.69;    % z-axis secular frequency, Hz

% calculation of a and q Mathieu parameters
% a_x = 2*wz^2/RF^2;   
% a0 = 2*(wx^2 - wy^2)/RF^2;
% qx = sqrt(4*(wx^2+wy^2+wz^2)/RF^2);
% qy = -qx;
% 
% ax = a0-a_x;
% ay = -a0-a_x;
% az = 2*a_x;

%% 
%   Create the LAMMPS experiment.
sim = LAMMPSSimulation();

% Add a simulation box. This determines the region that will be simulated.
% The simulation box may expand beyond these limits, but these specify the
% minimum dimensions of the simulation box.
sim.SetSimulationDomain(1e-2,1e-2,1e-2);
masses = [40 43 40 43 40 43 40 43 40 43 40];
% Add all atoms:
N = size(masses,2);  % number of ions
%type = [1 2 3 ]; % the type of ions
chars = ones(1, N);
% Set the index of ion with respect to which you want to calculate the modes.
ind = 1;
%trapmass = masses(ind);

wx = zeros(N, 1);   % column for x-axis secular frequencies
wy = zeros(N, 1);   % column for y-axis secular frequencies
wz = zeros(N, 1);   % column for axial secular frequencies


% Add ions to the simulation in the points with coordinates (0, 0, cor).
for i = 1:N
    ions = sim.AddAtomType(chars(i), masses(i));    % add i_th ion
    cor = (-N/2+i)*1e-5;    % z-coordinate of i_th ion
    corAtoms = placeAtoms(sim, ions, 0, 0, cor);    % place atoms to the cor
    % calculation of secular frequencies for i_th ion
    wx(i) = RF/2*sqrt(ax*masses(1)/masses(i)+(qx*masses(1)/masses(i))^2/2);
    wy(i) = RF/2*sqrt(ay*masses(1)/masses(i)+(qy*masses(1)/masses(i))^2/2);
    wz(i) = RF/2*sqrt(az*masses(1)/masses(i));
    sim.Add(my_PTrap(wx(i), wy(i), wz(i), ions));   % adds Paul trap in 
    % pseudopotential approximation, it's own trap for each ion
end

% Add a zero-temperature damping bath to reach equilibrium positions.
T = 0;
bath=langevinBath(T, 3e-5);
sim.Add(bath);

% Create files to write the results (coordinates, masses, charges)..
sim.Add(dump('positions.txt', {'id', 'mass', 'q', 'x', 'y', 'z'}, 1));
sim.Add(dump('secV.txt', {'id', timeAvg({'vx', 'vy', 'vz'}, 1/RF)}));

% Run simulation
sim.Add(evolve(20000));
sim.Execute();

%% Load and process the data
% Load the results from the output file:
[timesteps, id, mass, q, x,y,z] = readDump('positions.txt');
[~,~, sx, sy, sz] = readDump('secV.txt');
time = timesteps*sim.TimeStep;

% Check that the crystal remains in the same order.
for i=1:size(timesteps,2)
    if sumabs(round(mass(:,i)/amu)' - masses)>0
        error('The form of the crystal has changed (masses)! Fix it.')
    end

end
if sumabs(round(q(:,end)/ech)' - chars)>0
    error('The form of the crystal has changed (charges)! Fix it.')
end

% calculation of mass^(-0.5) matrix.
mss = zeros(1,N*3);
for i = 1:N
    mss(i:N:end) = [mass(i)^(-0.5) mass(i)^(-0.5) mass(i)^(-0.5)];
end
massmat = diag(mss);
l = (((chars(ind) * ech).^2)/(4*pi*eps0)/(masses(ind) * amu *((wz(ind)*2*pi)^2))).^(1/3);

% Get spectrum of each normal mode
u = x(:, end)/l;    
v = y(:, end)/l;
w = z(:, end)/l;
am = zeros(N*3);    % matrix for Hessian calculation
mw2 = masses(ind)*amu*(wz(ind)*2*pi)^2; 

for i=1:N
am(i,i)=mass(i)*(wx(i)*2*pi)^2/mw2;
         am(i+N,i+N)=mass(i)*(wy(i)*2*pi)^2/mw2;
         am(i+2*N,i+2*N)=mass(i)*(wz(i)*2*pi)^2/mw2;
    for j=1:N
       if i ~= j
        dij = sqrt((u(i)-u(j))^2+(v(i)-v(j))^2+(w(i)-w(j))^2)/chars(i)/chars(j);
        am(i,i) = am(i,i)-1/dij^3+3*(u(i)-u(j))^2/dij^5;
        am(i+N,i+N) = am(i+N,i+N)-1/dij^3+3*(v(i)-v(j))^2/dij^5;
        am(i+2*N,i+2*N) = am(i+2*N,i+2*N)-1/dij^3+3*(w(i)-w(j))^2/dij^5;
        am(i,i+N) = am(i,i+N)+3*(u(i) - u(j))*(v(i)-v(j))/dij^5;
        am(i+N, i) = am(i,i+N);
        am(i,i+2*N) = am(i,i+2*N)+3*(u(i) - u(j))*(w(i)-w(j))/dij^5;
        am(i+2*N, i) = am(i,i+2*N);
        am(i+N,i+2*N) = am(i+N,i+2*N)+3*(w(i) - w(j))*(v(i)-v(j))/dij^5;
        am(i+2*N, i+N) = am(i+N,i+2*N);
          am(i,j) = 1/dij^3-3*(u(i)-u(j))^2/dij^5;
          am(j,i) = am(i,j);
          am(i+N,j+N) = 1/dij^3-3*(v(i)-v(j))^2/dij^5;
          am(j+N, i+N) = am(i+N,j+N);
          am(i+2*N,j+2*N) = 1/dij^3-3*(w(i)-w(j))^2/dij^5;
          am(j+2*N,i+2*N) = am(i+2*N,j+2*N);
          am(i,j+N) = 3*(u(i) - u(j))*(v(j)-v(i))/dij^5;
          am(j+N, i) = am(i,j+N);
          am(i,j+2*N) = 3*(u(i) - u(j))*(w(j)-w(i))/dij^5;
          am(j+2*N, i) = am(i,j+2*N);
          am(i+N,j+2*N) = 3*(w(i) - w(j))*(v(j)-v(i))/dij^5;
          am(j+2*N,i+N) = am(i+N,j+2*N);
       end
    end
end
arm = am*masses(ind)*amu;
spmodes = massmat*arm*massmat;

[modes, spfreqs]=eig(spmodes);
     spfrsq=diag(spfreqs);
     freqs = spfrsq.^0.5;
     ansfrs = freqs
     ansmodes = modes;

plot(time, z)





     
     
     
     