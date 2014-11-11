function [mean_pop_firing_rate, mean_pop_norm_spike_pairs, Vs_traces, Vd_traces, s_traces] = gj_inhib_input(T0,no_cells,p_inhib,max_j,p_gj,gj_strength)
tic

%global dt t no_e_inputs e_rate e_size gj_strength no_i_inputs i_rate i_size CE_e CE_i max_k epsp ipsp rect
dt = .005; %this is in milliseconds
T = floor(T0/dt);
t = (1:T)*dt;

no_e_inputs = 127*no_cells; %127 = number of AMPA input synapses per cell in Hjorth et al, times 10 cells
e_rate = 2; % presynaptic firing rate (Hz) in Hjorth et al
e_size = 0.0053; %changes magnitude of input. 0.0053 gets you between 5 and 2 Hz firing rate.
%you stop getting firing rate decreases as conductance increases for 10 cells around 0.015. 
%at around 0.5 you can get increased spike pairs with increased conductance in 10 cell networks, but firing rate is deeply weird

inhib_strength = 3*(5*i_size)/sqrt(no_cells); %magnitude of steps of strength of gap junction

no_i_inputs = 93*no_cells; % = 93 times 10
i_rate = 2;
i_size = 0.0053;

if nargin < 6
    gj_strength = (5*e_size)/sqrt(no_cells);
end
if nargin < 5
    p_gj = 0;
end

tau_i1 = 1; tau_ir = 0.5; tau_id = 5; tau_i = 10; tau_r = 1;
tau_e1 = 1; tau_er = 0.5; tau_ed = 2;

CE_e = repmat(eye(no_cells), 1, no_e_inputs/no_cells);
CE_i = repmat(eye(no_cells), 1, no_i_inputs/no_cells);    %Define connectivity from inputs to cells.

max_k = 11; %number of conductance values to use

% EPSP for spikes at time t = 0.
epsp = tau_i*(exp(-max(t - tau_e1,0)/tau_ed) - exp(-max(t - tau_e1,0)/tau_er))/(tau_ed - tau_er);
epsp = epsp(epsp > eps);    %?
epsp = [zeros(1,length(epsp)) epsp]; %?

% IPSP for spikes at time t = 0.
ipsp = tau_i*(exp(-max(t - tau_i1,0)/tau_id) - exp(-max(t - tau_i1,0)/tau_ir))/(tau_id - tau_ir);
ipsp = ipsp(ipsp > eps);
ipsp = [zeros(1,length(ipsp)) ipsp];

%rectangle for convolutions
rect = 3*ones(1,5/dt);
rect = rect(rect > eps);
rect = [zeros(1,length(rect)) rect];

firing_rate = zeros(no_cells, max_k, max_j);
spike_pairs = zeros(max_k, max_j);
%firing_avg = zeros(max_k,1);
%pair_avg = zeros(max_k, max_j);

Vs_traces = zeros(max_j, max_k ,no_cells,T);
Vd_traces = zeros(max_j, max_k, no_cells,T);
s_traces = zeros(max_j, max_k, no_cells,T);

params = struct('dt', dt, 't', t, 'T', T, 'no_e_inputs', no_e_inputs, 'e_rate', e_rate, 'e_size', e_size, 'gj_strength', gj_strength, 'no_i_inputs', no_i_inputs, 'i_rate', i_rate, 'i_size', i_size, 'CE_e', CE_e, 'CE_i', CE_i, 'max_k', max_k, 'epsp', epsp, 'ipsp', ipsp, 'T0', T0, 'no_cells', no_cells, 'p_gj', p_gj, 'rect', rect, 'p_inhib', p_inhib, 'inhib_strength', inhib_strength);


parfor j = 1:max_j
	j
    [firing_rate(:,:,j), spike_pairs(:,j), Vs_traces(j,:,:,:), Vd_traces(j,:,:,:), s_traces(j,:,:,:)] = trial(params);
end
%spike_pairs = sum(spike_pairs,1);
%spike_pairs = spike_pairs./(2*firing_rate*(T0/1000)); %element-wise
%pair_avg = sum(spike_pairs,3)/max_j


pop_firing_rate = reshape(sum(firing_rate), max_k, max_j)*(1000/(T0*no_cells));
%pop_norm_spike_pairs = reshape(sum(sum(norm_spike_pairs)), max_k, max_j);

%spike_total = zeros(max_k, max_j);
%for a = 1:no_cells
%	for b = 1:no_cells
%		spike_total = spike_total + squeeze(norm_spike_pairs(a,b,:,:));
%	end
%end

%spike_total
spike_pairs
mean_pop_firing_rate = mean(pop_firing_rate, 2)
spike_pairs = spike_pairs./(2*pop_firing_rate*(T0/1000)); %element-wise
mean_pop_norm_spike_pairs = sum(spike_pairs,2)/max_j

%mean_pop_norm_spike_pairs = mean(spike_total./pop_firing_rate, 2)
%spike_pairs = spike_pairs./(2*firing_rate); %element-wise
%pair_avg = sum(spike_pairs,2)/max_j
toc

save('inhibsmalldata.mat','','-v7')
save('inhibbigdata.mat','Vs_traces','Vd_traces','s_traces','-v7.3')

inhib_intervals = (0:inhib_strength:(max_k-1)*inhib_strength);
figure
plot(inhib_intervals,mean_pop_firing_rate)
str = ['Average firing rate, ',num2str(no_cells), ' cells, uncorrelated input, ' num2str(max_j), ' trials'];
title(str)
xlabel('Inhibitory conductance')
ylabel('Firing rate (Hz)')
savefig('10_inhib_fr.fig')

figure
plot(inhib_intervals,mean_pop_norm_spike_pairs)
str = ['Spike pairs, ',num2str(no_cells), ' cells, uncorrelated input, ' num2str(max_j), ' trials'];
title(str)
xlabel('Inhibitory conductance')
ylabel('Proportion of paired spikes')
savefig('10_inhib_pairs.fig')
end

function [firing_rate_jslice, norm_spike_pairs_jslice, Vs_traces_jslice, Vd_traces_jslice, s_traces_jslice] = trial(params)
	firing_rate_jslice = zeros(params.no_cells, params.max_k, 1);
	%norm_spike_pairs_jslice = zeros(params.no_cells, params.no_cells, params.max_k, 1);
	norm_spike_pairs_jslice = zeros(params.max_k, 1);
	Vs_traces_jslice = zeros(1, params.max_k, params.no_cells, params.T);
	Vd_traces_jslice = zeros(1, params.max_k, params.no_cells, params.T);
    s_traces_jslice = zeros(1, params.max_k, params.no_cells, params.T);
	
    e_spikes = rand(params.no_e_inputs,length(params.t));
    e_spikes = e_spikes < params.e_rate*params.dt/1000;
    
    e_spike_arrivals = params.CE_e*e_spikes; % Calculating presynaptic spikes for each cell.

    epsps = nan(size(e_spike_arrivals)); % Calculating EPSP experienced by each cell.
    for c = 1:params.no_cells
      epsps(c,:) = params.e_size*conv(e_spike_arrivals(c,:),params.epsp,'same');
    end
    
    i_spikes = rand(params.no_i_inputs,length(params.t));
    i_spikes = i_spikes < params.i_rate*params.dt/1000;
   
    i_spike_arrivals = params.CE_i*i_spikes; % Calculating presynaptic spikes for each cell.

    ipsps = nan(size(i_spike_arrivals)); % Calculating IPSP experienced by each cell.

    for c = 1:params.no_cells
        ipsps(c,:) = params.i_size*conv(i_spike_arrivals(c,:),params.ipsp,'same');
    end

	
	CG = params.gj_strength*(rand(params.no_cells) < params.p_gj);
	CG = triu(CG);
	CG = CG + CG.';
    CG(logical(eye(size(CG)))) = 0;
    
    CS = params.inhib_strength*(rand(params.no_cells) < params.p_inhib);
	CS(logical(eye(size(CS)))) = 0;
    
    parfor k = 1:params.max_k
		k
		[firing_rate_jslice(:,k,1), norm_spike_pairs_jslice(k,1), Vs_traces_jslice(1,k,:,:), Vd_traces_jslice(1,k,:,:), s_traces_jslice(1,k,:,:)] = gj_value(params, k, epsps, ipsps, CG, CS);
    end
    %firing_avg = sum(firing_rate,2)/max_j
end

function [firing_rate_kslice, spike_pairs_kslice, Vs_traces_kslice, Vd_traces_kslice, s_traces_kslice] = gj_value(params, k, epsps, ipsps, CG, CS)
	firing_rate_kslice = zeros(params.no_cells, 1, 1);
	%spike_pairs_kslice = zeros(params.no_cells, params.no_cells, 1, 1);
	spike_pairs_kslice = zeros(1,1);
	Vs_traces_kslice = zeros(1, 1, params.no_cells, params.T);
	Vd_traces_kslice = zeros(1, 1, params.no_cells, params.T);
    s_traces_kslice = zeros(1, 1, params.no_cells, params.T);

    spike_indicator = zeros(params.no_cells,(params.T0/params.dt)-1);
	synch_indicator = zeros(params.no_cells, params.no_cells, (params.T0/params.dt)-2);
    %delta_t = 5; %number of milliseconds between two spikes to consider them "synchronous"
    %synch_interval = ones(delta_t/params.dt, 1);
    [Vs,Vd_traces_kslice(1,1,:,:),s_traces_kslice(1,1,:,:),~,~,~,~] = ing_w_dendritic_gap_jxn(params.no_cells, epsps-ipsps, params.T0, [], (k-1)*CS, CG);
	Vs_traces_kslice(1,1,:,:) = Vs;
	
	for a = 1:params.no_cells
        Vs_pos = Vs > 0;
        Vs_sign_change = diff(Vs_pos(a,:), [], 2);
        spike_indicator(a,:) = Vs_sign_change == 1;
        firing_rate_kslice(a, 1, 1) = sum(spike_indicator(a, :));
        %spike_indicator(a, :) = conv(spike_indicator(a, :), synch_interval, 'same') > 0;
    end
	
	%for a = 1:params.no_cells
    %    for b = 1:params.no_cells
    %        if CG(a,b) > 0
    %            synch_indicator = spike_indicator(a,:) & spike_indicator(b,:);
    %            synch_indicator = diff(synch_indicator) == 1;
    %            spike_pairs_kslice(a, b, 1, 1) = sum(synch_indicator); %spike_pairs(k,j) = spike_pairs(k,j) + sum(synch_indicator);
    %        end
    %    end
    %end
	
	%norm_spike_pairs_kslice(:, :, 1, 1) = diag(1./sqrt(firing_rate_kslice(:, 1, 1)))*spike_pairs_kslice(:, :, 1, 1)*diag(1./sqrt(firing_rate_kslice(:, 1, 1)));
	%for i = (1:length(params.no_cells))
	%	for j = (1:length(params.no_cells))
	%		norm_spike_pairs_kslice(i, j, 1, 1) = spike_pairs_kslice(i, j, 1, 1)/firing_rate_kslice(i,1,1);
	%	end
	%end
		
    wide_spikes = zeros(params.no_cells,length(params.t)-1);
    for c = 1:params.no_cells
		wide_spikes(c,:) = conv(spike_indicator(c,:),params.rect,'same');
	end
			
	for b = 1:params.no_cells
		for c = 1:params.no_cells
			if CG(b,c) > 0 && b ~= c 
				wide_sum = wide_spikes(b,:) + wide_spikes(c,:);
				foo = wide_sum > 5;
				wide_synch = diff(foo, [], 2);
				synch_indicator(b,c,:) = wide_synch == 1;
				spike_pairs_kslice(1,1) = spike_pairs_kslice(1,1) + sum(synch_indicator(b,c,:));
			end
		end
	end
end