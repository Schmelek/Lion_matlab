function bars_norm_modes(masses, norm_modes)

N = size(norm_modes, 1)/3;

Ca_40_ind = [];
AncillaIon_ind = [];

for j=1:N
    if masses(j) == 40
        Ca_40_ind(end+1) = j;

    else
        AncillaIon_ind(end+1) = j;
    end
end

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