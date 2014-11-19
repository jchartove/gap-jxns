function [inhib_gj_fr, inhib_gj_pairs] = inhib_gj(T0,no_cells,p_gj,max_j,p_inhib)
smallvars = 1;
max_i = 11;
max_k = 11;
n = (150*0.0053)/sqrt(no_cells);

%inhib conductance vs gj conductance vs firing rate, spike pairs. 
%inhib conductance = 0 to 3x highest gj conductance, in same steps as gj
inhib_gj_fr = zeros(max_i, max_k);
inhib_gj_pairs= zeros(max_i, max_k);
for i = (1:max_i)
    inhib_strength = (i-1) * n;
    [mean_pop_firing_rate, mean_pop_norm_spike_pairs,~,~,~] = gj_uncorr_input(T0,no_cells,p_gj,max_j,smallvars,p_inhib,inhib_strength);
    inhib_gj_fr(i,:) = mean_pop_firing_rate; %each row is one inhibitory value; each column is one gap junction value
    inhib_gj_pairs(i,:) = mean_pop_norm_spike_pairs;
end

str = ['inhib_gj_data', num2str(T0), '_', num2str(no_cells),'_',num2str(p_gj),'_',num2str(max_j),'_',num2str(p_inhib),'_',num2str(n),'.mat'];
save(str,'','-v7')

figure
imagesc([0 (50*0.0053)/sqrt(no_cells)],[0 (max_i-1)*n], inhib_gj_fr)
colorbar
str = ['Average firing rate, ',num2str(no_cells), ' cells, uncorrelated input, inhibition, ' num2str(max_j), ' trials'];
title(str)
xlabel('Gap junction conductance')
ylabel('Inhibitory conductance')
savefig('inhib_gj_fr.fig')

figure
imagesc([0 (50*0.0053)/sqrt(no_cells)], [0 (max_i-1)*n], inhib_gj_pairs)
colorbar
str = ['Spike pairs, ',num2str(no_cells), ' cells, uncorrelated input, inhibition, ' num2str(max_j), ' trials'];
title(str)
xlabel('Gap junction conductance')
ylabel('Inhibitory conductance')
savefig('inhib_gj_pairs.fig')
end