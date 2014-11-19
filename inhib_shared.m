function [inhib_shared_fr, inhib_shared_pairs] = inhib_shared(T0,no_cells,p_gj,max_j,p_inhib)
smallvars = 1;
max_i = 11;
max_k = 11;
n = (150*0.0053)/sqrt(no_cells);
gj_strength = (70*0.0053)/sqrt(no_cells);

%inhib conductance vs shared input vs firing rate, spike pairs
%gj conductance = whatever it is with correlated input tests
inhib_shared_fr = zeros(max_i, max_k);
inhib_shared_pairs= zeros(max_i, max_k);
for i = (1:max_i)
    inhib_strength = (i-1) * n;
    [mean_pop_firing_rate, mean_pop_norm_spike_pairs,~,~,~] = gj_corr_input(T0,no_cells,p_gj,max_j,smallvars,p_inhib,inhib_strength,gj_strength);
    inhib_shared_fr(i,:) = mean_pop_firing_rate;
    inhib_shared_pairs(i,:) = mean_pop_norm_spike_pairs;
end

str = ['inhib_shared_data', num2str(T0), '_', num2str(no_cells),'_',num2str(p_gj),'_',num2str(max_j),'_',num2str(p_inhib),'_',num2str(n),'.mat'];
save(str,'','-v7')

figure
imagesc([0 100], [0 (max_i-1)*n], inhib_shared_fr)
colorbar
str = ['Average firing rate, ',num2str(no_cells), ' cells, correlated input, inhibition, ' num2str(max_j), ' trials'];
title(str)
xlabel('Percent of input shared')
ylabel('Inhibitory conductance')
savefig('inhib_shared_fr.fig')

figure
imagesc([0 100], [0 (max_i-1)*n], inhib_shared_pairs)
colorbar
str = ['Spike pairs, ',num2str(no_cells), ' cells, correlated input, inhibition, ' num2str(max_j), ' trials'];
title(str)
xlabel('Percent of input shared')
ylabel('Inhibitory conductance')
savefig('inhib_shared_pairs.fig')
end