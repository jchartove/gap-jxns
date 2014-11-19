function [gj_shared_fr, gj_shared_pairs] = gj_shared(T0,no_cells,p_gj,max_j,p_inhib)
smallvars = 1;
max_i = 11;
max_k = 11;
n = (150*0.0053)/sqrt(no_cells);

%might be interesting: gj conductance vs shared input vs firing rate, spike pairs
gj_shared_fr = zeros(max_i, max_k);
gj_shared_pairs= zeros(max_i, max_k);

inhib_strength = n;
for i = (1:max_i)
    gj_strength = (i-1) * n/10;
    [mean_pop_firing_rate, mean_pop_norm_spike_pairs,~,~,~] = gj_corr_input(T0,no_cells,p_gj,max_j,smallvars,p_inhib,inhib_strength,gj_strength);
    gj_shared_fr(i,:) = mean_pop_firing_rate;
    gj_shared_pairs(i,:) = mean_pop_norm_spike_pairs;
end

str = ['gj_shared_data', num2str(T0), '_', num2str(no_cells),'_',num2str(p_gj),'_',num2str(max_j),'_',num2str(p_inhib),'_',num2str(n),'.mat'];
save(str,'','-v7')

figure
imagesc([0 100], [0 (max_i-1)*n/10], gj_shared_fr) 
colorbar
str = ['Average firing rate, ',num2str(no_cells), ' cells, correlated input, ' num2str(max_j), ' trials'];
title(str)
xlabel('Percent of input shared')
ylabel('Gap junction conductance')
savefig('gj_shared_fr.fig')

figure
imagesc([0 100], [0 (max_i-1)*n/10], gj_shared_pairs)
colorbar
str = ['Spike pairs, ',num2str(no_cells), ' cells, correlated input, ' num2str(max_j), ' trials'];
title(str)
xlabel('Percent of input shared')
ylabel('Gap junction conductance')
savefig('gj_shared_pairs.fig')
end