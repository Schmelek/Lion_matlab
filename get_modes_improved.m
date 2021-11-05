function [x_eq, y_eq, z_eq, norm_modes, frs, w_ind, l_out] = get_modes_improved(masses, chars, RF, ax, qx)

ech = 1.602176634e-19;  % electron charge, C
amu = 1.66053906660e-27;    % atomic mass unit, kg
eps0 = 8.8541878128e-12;    % vacuum electric permittivity

qy = -qx;
ay = ax;
az = -2*ax;
sim = LAMMPSSimulation();
sim.SetSimulationDomain(1,1,1);

% Add all atoms:
N = size(masses,2);  % number of ions
%type = [1 2 3 ]; % the type of ions

% Set the index of ion with respect to which you want to calculate the modes.
ind = 2;
%trapmass = masses(ind);

wx = zeros(N, 1);   % column for x-axis secular frequencies
wy = zeros(N, 1);   % column for y-axis secular frequencies
wz = zeros(N, 1);   % column for axial secular frequencies

Ion_type_1st = sim.AddAtomtype(chars(1), masses(1));
Ion_type = [Ion_type_1st];
pointer_array = [1];

for i = 1:(N-1)
    if ismember(masses(i+1), masses(1:i)) == True
       tmp = find(masses(1:i) == masses(i+1));
       pointer_array(end+1) = tmp(1);
    else
       Ion_type(end+1) = sim.AddAtom(chars(i+1), masses(i+));
       pointer_array(end+1) = i+1;
    end
end

for s = 1:N
    cor = (-N/2+s)*1e-4;    % z-coordinate of i_th ion
    placeAtoms(sim, Ion_type(pointer_array(s)), 0, 0, cor);    % place atoms to the cor
    % calculation of secular frequencies for i_th ion
    wx(s) = RF/2*sqrt(ax*masses(2)/masses(s)+(qx*masses(2)/masses(s))^2/2);
    wy(s) = RF/2*sqrt(ay*masses(2)/masses(s)+(qy*masses(2)/masses(s))^2/2);
    wz(s) = RF/2*sqrt(az*masses(2)/masses(s));
    sim.Add(my_PTrap(wx(s), wy(s), wz(s), Ion_type(pointer_array(s))));   % adds Paul trap in 
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
    frs = freqs;
    ansmodes = modes;
    norm_modes = ansmodes;
    
    tmp_norm_modes = [frs'; norm_modes];
    tmp_norm_modes = sortrows(tmp_norm_modes', 2+2*N, 'descend', 'ComparisonMethod','abs')';
    tmp_norm_modes = sortrows(tmp_norm_modes', 2+N, 'descend', 'ComparisonMethod','abs')';
    tmp_norm_modes = sortrows(tmp_norm_modes', 2, 'descend', 'ComparisonMethod','abs')';
    norm_modes = tmp_norm_modes(2:end, :);
    frs = tmp_norm_modes(1, :)'; 

    z_eq = (z(:, end));
    x_eq = vpa(x(:, end));
    y_eq = vpa(y(:, end));
    w_ind = wz(ind);
    l_out = l;
    if isreal(frs) == 0
        error('Frequencies are complex')
    end