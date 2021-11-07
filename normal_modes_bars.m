clearvars
sim.GPUAccel = 1;

% Physical constants

ech = 1.602176634e-19;  % electron charge, C
amu = 1.66053906660e-27;    % atomic mass unit, kg
eps0 = 8.8541878128e-12;    % vacuum electric permittivity

Vac = 850;  % AC-voltage, V
Udc = 3;   % DC radial voltage, V
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

masses = [43 43 43 40 40 40 40 40 40 40];
chars = ones(1, size(masses, 2));
N = size(masses, 2);

[x_eq, y_eq, z_eq, norm_modes, frs, w_n, l] = get_modes(masses, chars, RF, ax, qx);

%% Bars
for i=1:2*N
    subplot(N,2,i)
    if mod(i,2) == 1 
        if i == 1
            b = bar(norm_modes(1:N,(i+1)/2), 'FaceColor', 'flat');
            b.CData(Ca_40_ind, :) = repmat([1 0 0], size(Ca_40_ind, 2),1);
            b.CData(AncillaIon_ind, :) = repmat([0 0 1], size(AncillaIon_ind, 2),1);
            set(gca, 'YTick',-1:1:1,'Fontsize', 1);
            title('Radial(x) normal modes', 'Fontsize', 22)
        else
            b = bar(norm_modes(1:N,(i+1)/2), 'FaceColor', 'flat');
            b.CData(Ca_40_ind, :) = repmat([1 0 0], size(Ca_40_ind, 2),1);
            b.CData(AncillaIon_ind, :) = repmat([0 0 1], size(AncillaIon_ind, 2),1);
            set(gca, 'YTick',-1:1:1, 'Fontsize', 1);
        end
    else
        if i == 2
            b = bar(norm_modes((2*N+1):(3*N),2*N + i/2), 'FaceColor', 'flat');
            b.CData(Ca_40_ind, :) = repmat([1 0 0], size(Ca_40_ind, 2),1);
            b.CData(AncillaIon_ind, :) = repmat([0 0 1], size(AncillaIon_ind, 2),1);
            set(gca, 'YTick',-1:1:1,'Fontsize', 1);
            title('Axial normal modes', 'Fontsize', 22)
        else
            b = bar(norm_modes((2*N+1):(3*N),2*N + i/2), 'FaceColor', 'flat');
            b.CData(Ca_40_ind, :) = repmat([1 0 0], size(Ca_40_ind, 2),1);
            b.CData(AncillaIon_ind, :) = repmat([0 0 1], size(AncillaIon_ind, 2),1);
            set(gca, 'YTick',-1:1:1,'Fontsize', 1);
        end
    end 
end
%% Heatmap
tick = 2;
subplot(1,2,1);
tmp = sortrows([frs(1:N)';norm_modes(1:N, 1:N)]', 'ascend')';
tmp_radial = tmp(2:end, :);
h = heatmap(tmp_radial','CellLabelColor','none');
colormap(redblue(200));
title('Radial normal modes');
grid off;
CustomXLabels = string(1:N);
CustomXLabels(mod(1:N,tick) ~= 0) = " ";
h.XDisplayLabels = CustomXLabels;
CustomYLabels = string(1:N);
CustomYLabels(mod(1:N,tick) ~= 0) = " ";
h.YDisplayLabels = CustomYLabels;
h.XLabel = 'Ion number';
h.YLabel = 'Mode number';
% subplot(1,4,2);
% tmp = sortrows([frs(1:N)';norm_modes(1:N, 1:N)]', 'descend')';
% barh(tmp(1, :), 1);
subplot(1,2,2)
tmp = sortrows([frs(2*N+1:3*N)';norm_modes(2*N+1:end, 2*N+1:3*N)]', 'ascend')';
tmp_axial = tmp(2:end, :);
h = heatmap(tmp_axial','CellLabelColor','none');
colormap(redblueu(200));
CustomXLabels = string(1:N);
CustomXLabels(mod(1:N,tick) ~= 0) = " ";
h.XDisplayLabels = CustomXLabels;
CustomYLabels = string(1:N);
CustomYLabels(mod(1:N,tick) ~= 0) = " ";
h.YDisplayLabels = CustomYLabels;
grid off;
title('Axial normal modes');
h.XLabel = 'Ion number';
h.YLabel = 'Mode number';
% subplot(1,4,4);
% tmp = sortrows([frs(2*N+1:3*N)';norm_modes(2*N+1:end, 2*N+1:3*N)]', 'descend')';
% barh(tmp(1, :), 1)
%%
tmp = sortrows([frs(1:N)';norm_modes(1:N, 1:N)]', 'descend')';
tmp_frs_radial = tmp(1, :);
barh(tmp_frs_radial)