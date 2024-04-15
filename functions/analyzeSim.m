function results = analyzeSim(model,lbl, isLIF, VPM_strategy, VPM_onset)
    % get FR per layer
    cutTime=100;
    Sampling_period = 0.01; % in ms
    Sampling_rate = 1/(Sampling_period/1000); % in Hz
    cutTime_vec = cutTime:Sampling_period:model.time(end);
    binSize = 2; % in ms
    boundaries = [-300 300];
    
    field_names = cellfun(@(x) split(x,'_'), fieldnames(model.model.fixed_variables),'uniformoutput',false);
    
    if nargin<2
        lbl=[];
        isLIF = 0;
        VPM_strategy = 'threshold';
        VPM_onset = [];
    elseif nargin<3
        isLIF = 0;
        VPM_strategy = 'threshold';
        VPM_onset = [];
    elseif nargin<4
        VPM_strategy = 'threshold';
        VPM_onset = [];
    elseif nargin<5
        VPM_strategy = 'threshold';
        VPM_onset = [];
    end
    
    pops = ['VPM';'L4E';'L4I';'L3E';'L3I';'L2E';'L2I'];
    try
        model=dsCalcFR(model);
    catch ME
    end
    
    results = []; %table([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],'variablenames',...
%         {'pop','label','time','spike_times','netcons','power_welch','freq_welch','vm_autox','vm_autox_lag','periodicity_per_cell',...
%         'xcorr_with_vpm_per_cell','xcorr_with_vpm_per_cell_bins','ETAs','STAs','vm_distrib','vm_distrib_bins'});
    
    try
        vpm_vm = model.VPM_V((cutTime/Sampling_period)+1:end,:);
        vpm_vm_full = model.VPM_V;
    catch ME
        combinedVPMpops = [model.VPMcorr_V model.VPMuncorr_V];
%         combinedVPMpops=combinedVPMpops(:, randperm(size(combinedVPMpops, 2)));
        vpm_vm = combinedVPMpops((cutTime/Sampling_period)+1:end,:);
        vpm_vm_full = combinedVPMpops;
    end
    vpm_psth = plotPsthAndRaster('psth',getSpikeTimes(vpm_vm,cutTime_vec),binSize,cutTime_vec(end),0,'k');
    
    for pp = 1:length(pops)
        pop = pops(pp,:);
        
        % connectivity matrices
        netcons = [];
        for fnms = 1:length(field_names)
            if strcmp(field_names{fnms, 1}{1, 1}, pop) && strcmp(field_names{fnms, 1}{end, 1}, 'netcon')
                netcons.(cell2mat(join(field_names{fnms, 1},'_'))) = model.model.fixed_variables.(cell2mat(join(field_names{fnms, 1},'_')));
            end
        end
        
        if contains(pop,'VPM')
            try 
                spk_times = [model.VPMcorr_V_spike_times model.VPMuncorr_V_spike_times];
            catch ME
                spk_times = model.VPM_V_spike_times;
            end
            Vm = vpm_vm;
            Vm_full = vpm_vm_full;
        else
            if ~isLIF
                spk_times = model.([pop '_V_spike_times']);
            else
                spiking_threshold = model.model.parameters.([pop '_thresh']) -1;
                reset = model.model.parameters.([pop '_reset']) +1;
                spk_times = arrayfun(@(x) getLIFspikes(model.([pop '_V'])(:,x),model.time,cutTime,spiking_threshold,reset,0), 1:size(model.([pop '_V']),2),'uniformoutput',false);
            end
            Vm = model.([pop '_V'])((cutTime/Sampling_period)+1:end,:);
            Vm_full = model.([pop '_V']);
        end
        % isi per cell
%         isi = cellfun(@(x) diff(x), spk_times,'uniformoutput',false)';
        
        % spike rate per cell
%         fr = cellfun(@(x) numel(x)./((model.time(end)-cutTime)/1000), spk_times,'uniformoutput',false)';
        
        % corr and power
        pxx=[]; fq=[];
        if ~any(any(isnan(Vm))) && ~any(isnan(mean(Vm,1))) %&& ~any(mean(Vm,1) > -20)
            % power from single Vm traces
            [pxx, fq] = pwelch(Vm,[],[],[],Sampling_rate, 'power');
            fq = fq((fq >= 0) & (fq <= 1000),:); % take frequencies between 0 an 1 kHz
            pxx = pxx((fq >= 0) & (fq <= 1000),:); % shorten power vector to the size of freq vector

            % autocorrelations from single Vm PSTH
            vm_autox = []; vm_lags = [];
            psth_xcorr_single = []; psth_lags_single_xcorr = [];
            periodicity_per_cell=[];
            for neuron = 1:size(Vm,2)
                % autocorrelations of single Vm traces
                [vm_autox(neuron,:), vm_lags(neuron,:)] = xcorr(Vm(:,neuron)-mean(Vm(:,neuron)), 'coeff');
                
                % calculate periodicity
                peakThreshold = 0.2;
                peakDistance = 15; %in ms, for vm
                lower = find(vm_lags(neuron,:) == boundaries(1)/Sampling_period);
                upper = find(vm_lags(neuron,:) == boundaries(2)/Sampling_period);
                neuronSig = vm_autox(neuron,lower:upper); 
                try
                    periodicity_per_cell(neuron,1) = getPeriodicityFromAutox('vm', neuronSig, peakThreshold, peakDistance, Sampling_period);
                catch ME
                    periodicity_per_cell(neuron,1) = NaN;
                end
                % crosscorrelation between vpm psth and single neuron psth
                psth_single_neuron = plotPsthAndRaster('psth',spk_times(neuron),binSize,cutTime_vec(end),0,'k');
                [psth_xcorr_single(neuron,:), psth_lags_single_xcorr(neuron,:)] = xcorr(vpm_psth,psth_single_neuron, 'coeff');
            end
            % get only positive part of xcorr and autocorr
            vm_autox = vm_autox(:,round(size(vm_autox,2)/2):end);
            vm_lags = vm_lags(:,round(size(vm_lags,2)/2):end);
            psth_xcorr_single = psth_xcorr_single(:,round(size(psth_xcorr_single,2)/2):end);
            psth_lags_single_xcorr = psth_lags_single_xcorr(:,round(size(psth_lags_single_xcorr,2)/2):end);
            
            % find ETAs and STAs per cortical pop
            if ~strcmp(pop, 'VPM')
                minDistVPMevent = 20;
                binSize_events = 10;
                minVPMevent_height = 5*binSize_events*0.7;
                
                switch VPM_strategy
                    case 'threshold'
%                         VPM_info = {VPM_strategy, minVPMevent_height};
                        VPM_info = {'peak', minVPMevent_height};
                    case 'timestamp'
                        VPM_info = {VPM_strategy, VPM_onset};
                end
                
                if isLIF
                    thresh_pre = model.model.parameters.([pop '_thresh']) -1;
                    reset_pre = model.model.parameters.([pop '_reset']) +1;
                else
                    thresh_pre = NaN;
                    reset_pre = NaN;
                end
                ETA_analysis = getSTA_ETA('ETA',vpm_vm_full,Vm_full,[],binSize_events,VPM_info,minDistVPMevent,0,cutTime_vec(end),isLIF,[], ...
                        [thresh_pre, reset_pre]... % threshold and reset for presynaptic pop
                        );
                
                postSynPops = findPostSynPops(pop, model);
                STA_analysis = [];
                for postPop = 1:length(postSynPops)
                    if contains(pop, 'E')
                        con_mat = model.model.fixed_variables.([postSynPops{postPop} '_' pop '_iAMPA_netcon']);
                    elseif contains(pop,'I')
                        con_mat = model.model.fixed_variables.([postSynPops{postPop} '_' pop '_iGABAa_netcon']);
                    end
                    if isLIF
                        thresh_post = model.model.parameters.([postSynPops{postPop} '_thresh']) -1;
                        reset_post = model.model.parameters.([postSynPops{postPop} '_reset']) +1;
                    else
                        thresh_post = NaN;
                        reset_post = NaN;
                    end
                    
                    STA_analysis{postPop,1} = getSTA_ETA('STA',model.([pop '_V']), model.([postSynPops{postPop} '_V']),con_mat,[],[],[],0,cutTime_vec(end),isLIF, ...
                        [thresh_pre, reset_pre],... % threshold and reset for presynaptic pop
                        [thresh_post, reset_post]); % threshold and reset for postsynaptic pop
                    STA_analysis{postPop,2} = [pop '->' postSynPops{postPop}];
                end
            else
                % dummy results for VPM
                ETA_analysis = [];
                STA_analysis = [];
            end
%             figure;
%             plotPsthAndRaster('psth',getSpikeTimes(vpm_vm,cutTime_vec),binSize,cutTime_vec(end),1,'k');
%             hold on
%             xlim([100 600])
%             cellfun(@(x) plot(x, [15 15], '-r'), ev_times)

            % get Vm distribution
            [vm_dist, vm_dist_bins] = getVmDistribution(Vm);
            
            results = [table({pop},{lbl},{model.time},{spk_times'},{netcons},{pxx},{fq},...
                {[]},{[]},{[]},...{vm_autox},{vm_lags},{periodicity_per_cell},...
                {psth_xcorr_single},{psth_lags_single_xcorr},{ETA_analysis},...
                {[]},...{STA_analysis},
                {vm_dist},{vm_dist_bins},'variablenames',...
                {'pop','label','time','spike_times','netcons','power_welch','freq_welch','vm_autox','vm_autox_lag','periodicity_per_cell',...
                'xcorr_with_vpm_per_cell','xcorr_with_vpm_per_cell_bins','ETAs','STAs','vm_distrib','vm_distrib_bins'});results];
        else
            disp(['Flagged membrane potential trace. Did not analyze ' pop ' in ' lbl '. Check your data.'])
            reason = [];
            if any(any(isnan(Vm)))
                reason{end+1} = 'nans; ';
            elseif any(isnan(mean(Vm,1)))
                reason{end+1} = 'nans; ';
            elseif any(mean(Vm,1) > -20)
                reason{end+1} = 'depo.';
            end
            fprintf('Reason(s) for flagging: %s\n',cell2mat(reason));
%             pause
        end
    end
end

function postSynPops = findPostSynPops(preSynPop, model)
    connNames = fieldnames(model.model.fixed_variables);
    if contains(preSynPop,'E')
        connMxs = connNames(contains(connNames,[preSynPop '_iAMPA_netcon']));
        connMxs_split = split(connMxs,'_');
        postSynPops = connMxs_split(:,1)';
    elseif contains(preSynPop,'I')
        connMxs = connNames(contains(connNames,[preSynPop '_iGABAa_netcon']));
        connMxs_split = split(connMxs,'_');
        postSynPops = connMxs_split(:,1)';
    end
end
