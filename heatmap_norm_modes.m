function heatmap_norm_modes(masses, norm_modes, frs, w_n, x_eq, y_eq, z_eq) 
tick = 2;
N = size(norm_modes, 1)/3;
%% Equilibrium structure part

pointer_array = [1];

for i = 1:(N-1)
    if ismember(masses(i+1), masses(1:i)) == true
       pointer_array(end+1) = pointer_array(find(masses(i+1) == masses(1:i),1));
    else
       pointer_array(end+1) = max(pointer_array) + 1;
    end
end
tmp = max(pointer_array);
color_column_unit = (0:1/tmp:1)';
map_unit = [zeros(tmp, 2), color_column_unit(2:end)];
map = map_unit(pointer_array, :);

subplot('Position', [0.05 0.7 0.9 0.25])
scatter3(z_eq, x_eq, y_eq, 80,map, 'filled')
view(-45, 60);
grid on;
title('Equilibrium structure');
xlabel('z, m');
ylabel('x, m');
zlabel('y, m');
%% Normal modes part

subplot('Position', [0.05 0.1 0.25 0.5]);
tmp = sortrows([frs(1:N)';norm_modes(1:N, 1:N)]', 'ascend')';
tmp_radial = tmp(2:end, :);
h = heatmap(tmp_radial','CellLabelColor','none');
colormap(redblueu(200));
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

subplot('Position', [0.37 0.1 0.12 0.5]);
tmp = sortrows([frs(1:N)';norm_modes(1:N, 1:N)]', 'descend')';
b = barh(tmp(1, :)*w_n, 1);
set(gca, 'YTick', []);
ylabel('Mode number'); 
xlabel('Frequency, Hz');

subplot('Position', [0.525 0.1 0.25 0.5]);
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

subplot('Position', [0.855 0.1 0.12 0.5]);
tmp = sortrows([frs(2*N+1:3*N)';norm_modes(2*N+1:end, 2*N+1:3*N)]', 'descend')';
b = barh(tmp(1, :)*w_n, 1);
set(gca, 'YTick', []);
ylabel('Mode number'); 
xlabel('Frequency, Hz');

