function [inhib_shared_fr, inhib_shared_pairs] = inhib_shared(T0,no_cells,p_gj,max_j,p_inhib)
max_i = 10;
max_k = 11;
n = (15*0.0053)/sqrt(no_cells);

%inhib conductance vs shared input vs firing rate, spike pairs
%gj conductance = whatever it is with correlated input tests
inhib_shared_fr = zeros(max_i, max_k);
inhib_shared_pairs= zeros(max_i, max_k);
parfor i = (1:max_i)
    inhib_strength = i * n;
    [mean_pop_firing_rate, mean_pop_norm_spike_pairs] = gj_corr_input_small(T0,no_cells,p_gj,max_j,p_inhib,inhib_strength);
    inhib_shared_fr(i,:) = mean_pop_firing_rate;
    inhib_shared_pairs(i,:) = mean_pop_norm_spike_pairs;
end

save('inhib_shared_data.mat','','-v7')

figure
imagesc(inhib_shared_fr)
str = ['Average firing rate, ',num2str(no_cells), ' cells, correlated input, inhibition, ' num2str(max_j), ' trials'];
title(str)
xlabel('Percent of input shared')
ylabel('Inhibitory conductance')
savefig('inhib_shared_fr.fig')

figure
imagesc(inhib_shared_pairs)
str = ['Spike pairs, ',num2str(no_cells), ' cells, correlated input, inhibition, ' num2str(max_j), ' trials'];
title(str)
xlabel('Percent of input shared')
ylabel('Inhibitory conductance')
savefig('inhib_shared_pairs.fig')
end