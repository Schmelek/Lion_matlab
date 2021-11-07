function [x_eq, y_eq, z_eq, norm_modes, frs, w_ind, l_out] = get_modes(masses, chars, RF, ax, qx)
ech = 1.602176634e-19;  % electron charge, C
amu = 1.66053906660e-27;    % atomic mass unit, kg
eps0 = 8.8541878128e-12;    % vacuum electric permittivity

qy = -qx;
ay = ax;
az = -2*ax;

sim = LAMMPSSimulation();
sim.SetSimulationDomain(1e-2,1e-2,1e-2);

% Add all atoms:
N = size(masses,2);  % number of ions
%type = [1 2 3 ]; % the type of ions

% Set the index of ion with respect to which you want to calculate the modes.
ind = find(masses == 40, 1);
%trapmass = masses(ind);

%% Atoms and traps initializing
wx = zeros(N, 1);   % column for x-axis secular frequencies
wy = zeros(N, 1);   % column for y-axis secular frequencies
wz = zeros(N, 1);   % column for axial secular frequencies

Ion_type_1st = sim.AddAtomType(chars(1), masses(1));
Ion_type = [Ion_type_1st];
pointer_array = [1];

for i = 1:(N-1)
    if ismember(masses(i+1), masses(1:i)) == true
       pointer_array(end+1) = pointer_array(find(masses(i+1) == masses(1:i),1));
    else
       Ion_type(end+1) = sim.AddAtomType(chars(i+1), masses(i+1));
       pointer_array(end+1) = max(pointer_array) + 1;
    end
end

if mod(N,2) == 1
    cor = -(N-1)/2:1:(N-1)/2;
else
    cor = -N/2:1:N/2;
end

wx = RF/2*sqrt(ax*masses(ind)./masses+(qx*masses(ind)./masses).^2/2);
wy = RF/2*sqrt(ax*masses(ind)./masses+(qx*masses(ind)./masses).^2/2);
wz = RF/2*sqrt(az*masses(ind)./masses);

for i =1:N
    ions = placeAtoms(sim, Ion_type(pointer_array(i)), 0, 0, cor(i)*1e-5);
end

for i=1:size(Ion_type, 2)
    sim.Add(my_PTrap(wx(find(pointer_array == i,1)), wy(find(pointer_array == i,1)), wz(find(pointer_array == i,1)), Ion_type(i)));
end
    
%% Langevin bath & execution
T = 0;
bath=langevinBath(T, 3e-6);
sim.Add(bath);

sim.Add(dump('positions.txt', {'id', 'mass', 'q', 'x', 'y', 'z'}, 1));
sim.Add(dump('secV.txt', {'id', timeAvg({'vx', 'vy', 'vz'}, 1/RF)}));

% Run simulation
sim.Add(evolve(100000));
sim.Execute();
%% Normal modes
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
    x_eq = (x(:, end));
    y_eq = (y(:, end));
    w_ind = wz(ind);
    l_out = l;
    if isreal(frs) == 0
        error('Frequencies are complex')
    end
