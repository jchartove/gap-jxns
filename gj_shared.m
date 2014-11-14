function [gj_shared_fr, gj_shared_pairs] = gj_shared(T0,no_cells,p_gj,max_j,p_inhib)
max_i = 10;
max_k = 11;
n = (15*0.0053)/sqrt(no_cells);

%might be interesting: gj conductance vs shared input vs firing rate, spike pairs
gj_shared_fr = zeros(max_i, max_k);
gj_shared_pairs= zeros(max_i, max_k);

inhib_strength = n;
parfor i = (1:max_i)
    gj_strength = i * n;
    [mean_pop_firing_rate, mean_pop_norm_spike_pairs] = gj_corr_input_small(T0,no_cells,p_gj,max_j,p_inhib,inhib_strength,gj_strength);
    gj_shared_fr(i,:) = mean_pop_firing_rate;
    gj_shared_pairs(i,:) = mean_pop_norm_spike_pairs;
end

save('gj_shared_data.mat','','-v7')

figure
imagesc(gj_shared_fr) 
str = ['Average firing rate, ',num2str(no_cells), ' cells, correlated input, ' num2str(max_j), ' trials'];
title(str)
xlabel('Percent of input shared')
ylabel('Gap junction conductance')
savefig('gj_shared_fr.fig')

figure
imagesc(gj_shared_pairs)
str = ['Spike pairs, ',num2str(no_cells), ' cells, correlated input, ' num2str(max_j), ' trials'];
title(str)
xlabel('Percent of input shared')
ylabel('Gap junction conductance')
savefig('gj_shared_pairs.fig')
end