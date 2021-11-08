% clear the workspace

clearvars
sim.GPUAccel = 1;

% Physical constants

ech = 1.602176634e-19;  % electron charge, C
amu = 1.66053906660e-27;    % atomic mass unit, kg
eps0 = 8.8541878128e-12;    % vacuum electric permittivity



RF = 6e6;
qx = 0.06;  
ax = -0.0007;  
qy = qx;
ay = ax;
az = -2*ax;


minimisationSteps = 10000;
interval = 3000000;

array = [[40 40 40 40 40 43 43];[43 40 40 40 40 40 43]]
%[[40 40 40 23 40];[40 40 40 25 40]; [40 40 40 30 40];[40 40 40 33 40];[40 40 40 35 40];[40 40 40 36 40];[40 40 40 37 40] ;[40 40 40 38 40];[40 40 40 39 40]; [40 40 40 41 40];[40 40 40 42 40];[40 40 40 43 40];[40 40 40 44 40];[40 40 40 45 40];[40 40 40 50 40]; [40 40 40 60 40]; [40 40 40 80 40];[40 40 40 120 40]]; 


Array_size = size(array,1); 

% T_Ca = ones(Array_size, 28801);
% T_AncillaIon = ones(Array_size, 28801);
% T_Ca_x = ones(Array_size, 28801);
% T_Ca_y = ones(Array_size, 28801);
% T_Ca_z = ones(Array_size, 28801);
% T_Ca_transverse = ones(Array_size, 28801);
% 
% z_eq = ones(60); %number of permutations of every ancilla ion type * number of ancilla ion types * number of ions in a string
% norm_modes = ones(180, 15); %(3 * number of ions in a string * number of ancilla ion types * number of permutations),(3 * number of ions in a string)
% frs = ones(180, 1); %(3 * number of ions in a string * number of ancilla ion types * number of permutations)
% w_n = ones(12, 1); %number of ancilla ion types * number of permutations
% l = ones(12, 1); %%number of ancilla ion types * number of permutations

T_Ca_x = [];
T_Ca_y = [];
T_Ca_z = [];
T_Ca_transverse = [];
T_Ca = [];
T_AncillaIon = [];
T_sole_ion_x = [];
T_sole_ion_z = [];
T_sole_ion_y = [];
T_sole_ion_transverse = [];
T_sole_ion = [];

%if i_norm_mode_index exceeds N, the displacement is 0
x_norm_mode_ind = [1];
y_norm_mode_ind = [6];
z_norm_mode_ind = [1];

z_eq = [];
x_eq = [];
y_eq = [];
norm_modes = [];
frs = [];
w_n = [];
l = [];
tic  

for i=1:size(array,1)
    tmp_x_eq=[];
    tmp_y_eq=[];
    tmp_z_eq=[];
    tmp_norm_modes = [];
    tmp_frs = [];
    tmp_w_n = [];
    tmp_l = [];

    Ca_40_ind = [];
    AncillaIon_ind = [];
    masses = array(i, :);
    chars = ones(1, size(masses, 2));
    N = size(masses, 2);


    wx = zeros(N, 1);   % column for x-axis secular frequencies
    wy = zeros(N, 1);   % column for y-axis secular frequencies
    wz = zeros(N, 1);   % column for axial secular frequencies

    [~,~,~,tmp_x_eq, tmp_y_eq, tmp_z_eq, tmp_norm_modes, tmp_frs, tmp_w_n, tmp_l] = get_modes(masses, chars, RF, ax, qx);
    x_eq = [x_eq; tmp_x_eq];
    y_eq = [y_eq; tmp_y_eq];
    z_eq = [z_eq; tmp_z_eq];
    norm_modes = [norm_modes; tmp_norm_modes];
    frs = [frs;tmp_frs];
    w_n = [w_n, tmp_w_n];
    l = [l; tmp_l];
    
    for mode_ind = 1:size(z_norm_mode_ind, 2) 

        
        sim = LAMMPSSimulation();
        sim.SetSimulationDomain(1e-2,1e-2,1e-2);

        for j=1:N
            if masses(j) == 40
                Ca40 = AddAtomType(sim, 1, 40);
                Ca_40_ind(end+1) = j;

            else
                AncillaIon = AddAtomType(sim, 1, masses(j));
                AncillaIon_ind(end+1) = j;
            end
        end


        for s = 1:N
            ions = sim.AddAtomType(chars(s), masses(s));    % add i_th ion
%             corAtoms = placeAtoms(sim, ions, tmp_norm_modes(s ,x_norm_mode_ind(mode_ind))*tmp_l*0.035 , 0*tmp_norm_modes(s+N ,N+y_norm_mode_ind(1))*tmp_l*0.001 , tmp_z_eq(s) + tmp_norm_modes(s + 2*N ,2*N+z_norm_mode_ind(mode_ind))*tmp_l*0.035);
            % calculation of secular frequencies for i_th ion
            corAtoms = placeAtoms(sim, ions, 0, 0, (-N/2+s)*1e-5);
            wx(s) = RF/2*sqrt(ax*masses(2)/masses(s)+(qx*masses(2)/masses(s))^2/2);
            wy(s) = RF/2*sqrt(ay*masses(2)/masses(s)+(qy*masses(2)/masses(s))^2/2);
            wz(s) = RF/2*sqrt(az*masses(2)/masses(s));
            sim.Add(my_PTrap(wx(s), wy(s), wz(s), ions));   % adds Paul trap in 
            % pseudopotential approximation, it's own trap for each ion
        end

        allBath = langevinBath(0.1, 3e-6); 
        sim.Add(allBath);
       
        sim.Add(evolve(minimisationSteps));
        sim.Remove(allBath);
        sim.Add(langevinBath(1e-3, 5e-7, sim.Group(AncillaIon_ind)));
        sim.Add(dump('sympcool.txt', {'id','mass' , 'x', 'y', 'z','vx', 'vy', 'vz'}, 20));

        sim.Add(evolve(interval));

        sim.Execute();

        [t, id,mass, x,y,z, sx,sy,sz] = readDump('sympcool.txt');
        t = (t-minimisationSteps)*sim.TimeStep;
        
        T_sole_ion_x = [T_sole_ion_x;sx(Ca_40_ind, :).^2* Const.amu * Ca40.Mass/3/Const.kB];
        T_sole_ion_y = [T_sole_ion_y;sy(Ca_40_ind, :).^2* Const.amu * Ca40.Mass/3/Const.kB];
        T_sole_ion_z = [T_sole_ion_z;sz(Ca_40_ind, :).^2* Const.amu * Ca40.Mass/3/Const.kB];
        T_Ca_x = [T_Ca_x;sum(sx(Ca_40_ind, :).^2,1)* Const.amu * Ca40.Mass/3/Const.kB/size(Ca_40_ind,2)];%/(tmp_w_n*tmp_frs(x_norm_mode_ind(mode_ind)))^2];
        T_Ca_y = [T_Ca_y;sum(sy(Ca_40_ind, :).^2,1)* Const.amu * Ca40.Mass/3/Const.kB/size(Ca_40_ind,2)];%/(tmp_w_n*tmp_frs(N+y_norm_mode_ind(mode_ind)))^2];
        T_Ca_z = [T_Ca_z;sum(sz(Ca_40_ind, :).^2,1)* Const.amu * Ca40.Mass/3/Const.kB/size(Ca_40_ind,2)];%/(tmp_w_n*tmp_frs(2*N+z_norm_mode_ind(mode_ind)))^2];
  
        T_AncillaIon = [T_AncillaIon; sum(sx(AncillaIon_ind, :).^2 + sy(AncillaIon_ind, :).^2 + sz(AncillaIon_ind, :).^2, 1)*AncillaIon.Mass*Const.amu/3/Const.kB/size(AncillaIon_ind,2)];

    end

end
        T_Ca = T_Ca_x + T_Ca_y + T_Ca_z;
        T_Ca_transverse = T_Ca_x + T_Ca_y;
        T_sole_ion = T_sole_ion_x + T_sole_ion_y +T_sole_ion_z;
        T_sole_ion_transverse = T_sole_ion_x + T_sole_ion_y;
toc
plot(t,z)
%% 
plot(t, T_Ca, 'b.')
grid on;
xlabel('t, s');
ylabel('Average temperature of Ca over the crystal, K');
title('40 40 40 40 40 40 41')
%% 
subplot(1,2,1)
plot(t, smoothdata(1e3*T_sole_ion_x(5,:), 'movmedian',5000), 'b',t, smoothdata(1e3*T_sole_ion_x(4,:), 'movmedian',5000), 'c', t, smoothdata(1e3*T_sole_ion_x(1,:), 'movmedian',5000),'r', 'Linewidth', 1.5)
grid on;
xlabel('t, s', 'Fontsize', 18);
ylabel('T, mK', 'Fontsize',18);
legend('5th ion','4th ion' ,'1st ion', 'Fontsize', 22);
title('Radial(x) temperatures of sole ions', 'Fontsize', 18);
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [0.5,1,0]);
set(gca, 'Fontsize', 20);
subplot(1,2,2)
plot(t(500:end), smoothdata(1e3*T_sole_ion_z(5,500:end), 'movmedian',5000), 'b',t(500:end), smoothdata(1e3*T_sole_ion_z(4,500:end), 'movmedian',5000), 'c', t(500:end), smoothdata(1e3*T_sole_ion_z(1,500:end), 'movmedian',5000),'r', 'Linewidth', 1.5)
grid on;
xlabel('t, s', 'Fontsize', 18);
ylabel('T, mK', 'Fontsize',18);
legend('5th ion','4th ion' ,'1st ion', 'Fontsize', 22);
title('Axial temperatures of sole ions', 'Fontsize', 18);
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [0.55,1,0]);
set(gca, 'Fontsize', 20);
%% 

for i=1:2*N
    subplot(N,2,i)
    if mod(i,2) == 1 
        if i == 1
            b = bar(norm_modes(3*N+1:4*N,(i+1)/2));
            b.FaceColor = 'flat';

            b.CData(2,:) = [10 0 0];
            b.CData(3,:) = [10 0 0];
            b.CData(4,:) = [10 0 0];
            b.CData(5,:) = [10 0 0];
            b.CData(6,:) = [10 0 0];
            set(gca, 'YTick',-1:1:1,'Fontsize', 18)
            title('Radial(x) normal modes', 'Fontsize', 22)
        else
                        b = bar(norm_modes(3*N+1:4*N,(i+1)/2));
            b.FaceColor = 'flat';

            b.CData(2,:) = [10 0 0];
            b.CData(3,:) = [10 0 0];
            b.CData(4,:) = [10 0 0];
            b.CData(5,:) = [10 0 0];
            b.CData(6,:) = [10 0 0];
            set(gca, 'YTick',-1:1:1, 'Fontsize', 18)
        end
    else
        if i == 2
            b = bar(norm_modes((5*N+1):(6*N),2*N + i/2));
            b.FaceColor = 'flat';

            b.CData(2,:) = [10 0 0];
            b.CData(3,:) = [10 0 0];
            b.CData(4,:) = [10 0 0];
            b.CData(5,:) = [10 0 0];
            b.CData(6,:) = [10 0 0];
            set(gca, 'YTick',-1:1:1,'Fontsize', 18)
            title('Axial normal modes', 'Fontsize', 22)
        else
            b = bar(norm_modes((5*N+1):(6*N),2*N + i/2));
            b.FaceColor = 'flat';

            b.CData(2,:) = [10 0 0];
            b.CData(3,:) = [10 0 0];
            b.CData(4,:) = [10 0 0];
            b.CData(5,:) = [10 0 0];
            b.CData(6,:) = [10 0 0];
            set(gca, 'YTick',-1:1:1,'Fontsize', 18)
        end
    end 
end
%% 
plot(t,smoothdata(T_sole_ion_x(6,:), 'movmedian',1000),'g', t,  T_sole_ion_x(6,:))
%% 
plot(t,z)
xlabel('t, s');
ylabel('z, m');
title('z(t) for 40 40 40 40 40 40 41');
%% 
plot(t, T_Ca)
grid on;

title('Sympathetic cooling of 40 40 40 40 41 crystal');
xlabel('time, s');
ylabel('T, K')
%% 
plot(t(2000:end), smoothdata(T_sole_ion_z(6,2000:end), 'movmedian',500))
%% 
[fourier, f] = make_F(t, x(1,:));

%% Data preparation for fiiting app
for i=1:(size(z_norm_mode_ind,2)*size(array,1))
    assignin('base', ['T_Ca_' num2str(i)], T_Ca(i, :));
    assignin('base', ['T_Ca_transverse_' num2str(i)], T_Ca_transverse(i, :));
    assignin('base', ['T_Ca_z_' num2str(i)], T_Ca_z(i, :));
end

for i=1:(N-1)
    assignin('base',['T_sole_ion_z' num2str(i)], T_sole_ion_z(i,:));
    assignin('base',['T_sole_ion_x' num2str(i)], T_sole_ion_x(i,:));
end;
tau = t(:);
%% Sequence of cooling rates for ion chains for large temperatures
Cooling_rate_side_5 = [2.468 3.145 4.44 6.452 7.298 7.501 7.905 8.338 8.771 9.254 9.309 9.365 9.286 8.123 8.469 7.471 2.116 0.4001];
CR_5_side = [1.608 2.209 3.987 4.804 5.891 6.457 6.919 7.178 7.187 7.148 7.256 7.277 7.263 7.098 5.781 3.5 2.532 0.5608]
CR_5_mid = [5.533 8.678 20.82 17.53 19.65 19.3 17.79 16.37 14.78 12.35 11.28 10.62 10.62 10.58 0.3869 0.4468 2.818 0.9707 ]
amuai = [23 25 30 33 35 36 37 38 39 41 42 43 44 45 50 60 80 120];
plot(amuai, CR_5_side, 'bo:')
%% Sequence of cooling rates for sole normal modes
CR_nm_t_41 = [0.2905 0.7819 0.9163 0.3157 0.7333];
CR_nm_a_41 = [0.3616 2.893 4.084 3.912 3.691];
CR_nm_a_41 = [0.4041 3.198 4.272 3.922 3.769];
CR_nm_t_41 = [0.4163 1.13 1.305 0.3963 1.183];
plot(frs(11:15), CR_nm_a_41, 'bo')
grid on;
hold on;
plot(frs(6:10), CR_nm_t_41, 'ro')
title('Cooling rates for sole normal modes for the crystal 40 40 40 40 41')
xlabel('Frequency, 1/s');
ylabel('Cooling rate, 1/s');
legend('Axial modes', 'Radial modes')
hold off;
%% Sequence of cooling rates for ion chains for small temperatures

Cooling_rate_side_5 = [29.5 80.34 169.4 230.5 267.6 297.3 308.6 317.1 318.8 312.8 306 298.4 290.5 250.8 189.5 125.4 76.6];
Cooling_rate_mid_5 = [97.65 291 629.7 823.6 909.4 953.3 960.3 960.4 926.2 895.8 864.2 832.3 800.7 660.7 478.4 312.7 194.5];
AncillaIon_amu = [20 25 30 33 35 37 38 39 41 42 43 44 45 50 60 80 120];
Cooling_rate_breathing_mode = zeros(17,1);
plot(AncillaIon_amu, Cooling_rate_side_5, 'bo:')
grid on;
hold on;
plot(AncillaIon_amu,Cooling_rate_mid_5, 'ro:');
hold on;
plot(AncillaIon_amu, Cooling_rate_breathing_mode, 'ko:')
xlabel('Ancilla Ion amu');
ylabel('Cooling rate, 1/s');
title('Sympathetic cooling for mixed species ion chains');
legend('40 40 40 40 AI', '40 40 40 AI 40', '40 40 AI 40 40');
hold off;
%%
Cooling_rate_side_6 = [40.33 78.76  141.9 182.7 207.4 228.6 237.5 245.1 248.9 245.4 241.3 236.5 231.4 204.4 160.3 111.7 72.79];
Cooling_rate_side_7 = [41.5 78.26 117.8 133.9 142.6 150.7 154.6 158.4 162 161.7 161.2 160.4 159.4 152 132.8 101.2 68.45];
Cooling_rate_1st_mode = [10000 14000 16000 15000 13500 10000];
AncillaIon_amu_mode_graph = [20 30 39 41 50 60]
plot(AncillaIon_amu, Cooling_rate_side_5, 'bo:')
hold on
plot (AncillaIon_amu, Cooling_rate_side_6, 'mo:')
hold on
plot (AncillaIon_amu, Cooling_rate_side_7, 'ro:')
hold on
xlabel('Ancilla Ion amu')
ylabel('Cooling rate, 1/s')
grid on;
title('Sympathetic cooling of mixed species ion chains')
legend('40 40 40 40 AI', '40 40 40 40 40 AI', '40 40 40 40 40 40 AI')
hold off
%% Calculating the minimal difference between normal mode frequencies
axial_frs = [];
transverse_frs = [];
for i=1:Array_size
    axial_frs = [axial_frs; frs(i*N*3-N+1:N*3*i)']; %the coefficients are unique for every permutation sequence
    transverse_frs = [transverse_frs; frs((i-1)*3*N+1:3*N*(i-1)+N)'];%the coefficients are unique for every permutation sequence
end

min_axial_frs_diff = min(abs([zeros(Array_size,1),sort(axial_frs,2)] - [sort(axial_frs,2), zeros(Array_size,1)]),[],2);
min_transverse_frs_diff = min(abs([zeros(Array_size,1),sort(transverse_frs,2)] - [sort(transverse_frs,2), zeros(Array_size,1)]),[],2);

plot(array(:, N), min_axial_frs_diff*wz(1), 'bo:')
hold on
plot(array(:, N), min_transverse_frs_diff*wz(1), 'ro:')
hold off
legend('axial', 'transverse')
xlabel('Ancilla Ion amu')
ylabel('minimal frequency difference, 1/s')