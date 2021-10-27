function [x_eq, y_eq, z_eq, norm_modes, ansfrs, w_ind, l_out] = get_modes(masses, chars, RF, ax, qx)
ech = 1.602176634e-19;  % electron charge, C
amu = 1.66053906660e-27;    % atomic mass unit, kg
eps0 = 8.8541878128e-12;    % vacuum electric permittivity

qy = qx;
ay = ax;
az = -2*ax;
sim = LAMMPSSimulation();
sim.SetSimulationDomain(1e-2,1e-2,1e-2);

% Add all atoms:
N = size(masses,2);  % number of ions
%type = [1 2 3 ]; % the type of ions

% Set the index of ion with respect to which you want to calculate the modes.
ind = 2;
%trapmass = masses(ind);

wx = zeros(N, 1);   % column for x-axis secular frequencies
wy = zeros(N, 1);   % column for y-axis secular frequencies
wz = zeros(N, 1);   % column for axial secular frequencies

for i = 1:N
    ions = sim.AddAtomType(chars(i), masses(i));    % add i_th ion
    cor = (-N/2+i)*1e-5;    % z-coordinate of i_th ion
    corAtoms = placeAtoms(sim, ions, 0, 0, cor);    % place atoms to the cor
    % calculation of secular frequencies for i_th ion
    wx(i) = RF/2*sqrt(ax*masses(2)/masses(i)+(qx*masses(2)/masses(i))^2/2);
    wy(i) = RF/2*sqrt(ay*masses(2)/masses(i)+(qy*masses(2)/masses(i))^2/2);
    wz(i) = RF/2*sqrt(az*masses(2)/masses(i));
    sim.Add(my_PTrap(wx(i), wy(i), wz(i), ions));   % adds Paul trap in 
    % pseudopotential approximation, it's own trap for each ion
end

T = 0;
bath=langevinBath(T, 3e-6);
sim.Add(bath);

sim.Add(dump('positions.txt', {'id', 'mass', 'q', 'x', 'y', 'z'}, 1));
sim.Add(dump('secV.txt', {'id', timeAvg({'vx', 'vy', 'vz'}, 1/RF)}));

% Run simulation
sim.Add(evolve(50000));
sim.Execute();

[timesteps, id, mass, q, x,y,z] = readDump('positions.txt');
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
    mss(i:N:end) = vpa([mass(i)^(-0.5) mass(i)^(-0.5) mass(i)^(-0.5)]);
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
     ansfrs = freqs;
     ansmodes = modes;
     norm_modes = vpa(ansmodes);
     z_eq = vpa(z(:, end));
     x_eq = vpa(x(:, end));
     y_eq = vpa(y(:, end));
     w_ind = wz(ind);
     l_out = l;