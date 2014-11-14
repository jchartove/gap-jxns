function [inhib_gj_fr, inhib_gj_pairs] = inhib_gj(T0,no_cells,p_gj,max_j,p_inhib)
max_i = 10;
max_k = 11;
n = (15*0.0053)/sqrt(no_cells);

%inhib conductance vs gj conductance vs firing rate, spike pairs. 
%inhib conductance = 0 to 3x highest gj conductance, in same steps as gj
inhib_gj_fr = zeros(max_i, max_k);
inhib_gj_pairs= zeros(max_i, max_k);
parfor i = (1:max_i)
    inhib_strength = i * n;
    [mean_pop_firing_rate, mean_pop_norm_spike_pairs] = gj_uncorr_input_small(T0,no_cells,p_gj,max_j,p_inhib,inhib_strength);
    inhib_gj_fr(i,:) = mean_pop_firing_rate;
    inhib_gj_pairs(i,:) = mean_pop_norm_spike_pairs;
end

save('inhib_gj_data.mat','','-v7')

figure
imagesc(inhib_gj_fr)
str = ['Average firing rate, ',num2str(no_cells), ' cells, uncorrelated input, inhibition, ' num2str(max_j), ' trials'];
title(str)
xlabel('Gap junction conductance')
ylabel('Inhibitory conductance')
savefig('inhib_gj_fr.fig')

figure
imagesc(inhib_gj_pairs)
str = ['Spike pairs, ',num2str(no_cells), ' cells, uncorrelated input, inhibition, ' num2str(max_j), ' trials'];
title(str)
xlabel('Gap junction conductance')
ylabel('Inhibitory conductance')
savefig('inhib_gj_pairs.fig')
end