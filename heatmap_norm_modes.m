function heatmap_norm_modes(norm_modes, frs) 
tick = 2;
N = size(norm_modes, 1)/3;

subplot('Position', [0.05 0.1 0.25 0.9]);
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

subplot('Position', [0.37 0.1 0.12 0.9]);
tmp = sortrows([frs(1:N)';norm_modes(1:N, 1:N)]', 'descend')';
b = barh(tmp(1, :), 1);
set(gca, 'YTick', []);

subplot('Position', [0.525 0.1 0.25 0.9]);
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

subplot('Position', [0.855 0.1 0.12 0.9]);
tmp = sortrows([frs(2*N+1:3*N)';norm_modes(2*N+1:end, 2*N+1:3*N)]', 'descend')';
b = barh(tmp(1, :), 1);
set(gca, 'YTick', []);

