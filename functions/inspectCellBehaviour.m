function cb = inspectCellBehaviour(model, cellToInspect, plon, sim_tag, threshold_for_inspection)
warning('off')
    try
        model=dsCalcFR(model);
    catch ME
    end
    time = [model.time(1) model.time(end) diff(model.time([1 2]))];        % time in ms in format [start stop step]
    cutStartTime =  100;  % in ms
    pops = ['L2E';'L2I';'L3E';'L3I';'L4E';'L4I';'VPM'];
    
    if ~any(string(pops)==cellToInspect)
        fprintf('%s is not a valid cell type\n',cellToInspect)
        cb=[];
        return
    end
    if nargin<5
        threshold_for_inspection = -55; % mV
    end
    
    % check for unsettled cell activity
    cellBehaviour = [];
    
    % remove all datapoints that overshoot (i.e., are above 0).. these
    % presumably belong to a spike, so i essentially remove
    % overshooting part of the spike
    faulty = find(cell2mat(arrayfun(@(x) mean(model.([cellToInspect '_V'])(model.([cellToInspect '_V'])(:,x)<0,x))>=threshold_for_inspection, 1:size(model.([cellToInspect '_V']),2),'uniformoutput',false)));
    healthy = find(cell2mat(arrayfun(@(x) mean(model.([cellToInspect '_V'])(model.([cellToInspect '_V'])(:,x)<0,x))<threshold_for_inspection, 1:size(model.([cellToInspect '_V']),2),'uniformoutput',false)));

    % see if there are spiking faulty cells
    faulty_not_spiking=[]; faulty_spiking=[];
    if ~isempty(faulty)
        for ii=1:length(faulty) 
            spikes_or_not = any(model.([cellToInspect '_V_spike_times']){1, faulty(ii)}>cutStartTime);
            if ~spikes_or_not
                faulty_not_spiking(end+1) = faulty(ii);
            else
                faulty_spiking(end+1) = faulty(ii);
            end
        end          
    end; clear ii
    cellBehaviour.(cellToInspect).fns = faulty_not_spiking;
    cellBehaviour.(cellToInspect).fs = faulty_spiking;
    cellBehaviour.(cellToInspect).ht = healthy;
   
    % plot faulty cells if the user wants to see
    if plon(1)==1 
        if isempty(faulty)
            fprintf('No faulty cells to plot.\n')
        else
            f1000=figure; tiledlayout(3,ceil(length(faulty)/3))
            for ii=1:length(faulty)
                if ismember(faulty(ii),faulty_not_spiking)
                    nexttile;
                    col='k';
                    plot(time(1):time(3):time(2),model.([cellToInspect '_V'])(:,faulty(ii)), 'color', col);
                    hold on;
                    plot([cutStartTime cutStartTime], [-100 50], ':k')
                    xlim([0 time(2)])
                    ylim([-100 50])
                    legend([cellToInspect ' cell idx: ' num2str(faulty(ii))])
                elseif ismember(faulty(ii),faulty_spiking)
                    nexttile;
                    col='b';
                    plot(time(1):time(3):time(2),model.([cellToInspect '_V'])(:,faulty(ii)), 'color', col);
                    hold on;
                    plot([cutStartTime cutStartTime], [-100 50], ':k')
                    xlim([0 time(2)])
                    ylim([-100 50])
                    legend([cellToInspect ' cell idx: ' num2str(faulty(ii))])
                end
            end; clear ii
        end
    end
    
    % find the names of presynaptic partners to each cell
    if contains(cellToInspect,'L2')
        presynPops = ["L2E", "L2I","L3E","L3I"];
    elseif contains(cellToInspect,'L3E')
        presynPops = ["L3E","L3I","L4E", "L4I"];
    elseif contains(cellToInspect,'L3I')
        presynPops = ["L3E","L3I","L4E"];
    elseif contains(cellToInspect,'L4')
        presynPops = ["L4E","L4I","VPM"];
    end
    % loop through presynaptic partners and count connections per cell
    PresynapticPartnerCount = [];
    for ii=1:length(presynPops)
        if contains(presynPops(ii),'I') 
            currentType = 'iGABAa';
        else
            currentType = 'iAMPA';
        end
        if contains(cellToInspect,'L4') && ii==3
            presynPartners = arrayfun(@(x) numel(find(model.model.fixed_variables.([char(cellToInspect) '_' char(presynPops(ii)) 'corr_' currentType '_netcon'])(:,x))) + ...
                numel(find(model.model.fixed_variables.([char(cellToInspect) '_' char(presynPops(ii)) 'uncorr_' currentType '_netcon'])(:,x))),...
                1:size(model.([cellToInspect '_V']),2),'uniformoutput',false);
            PresynapticPartnerCount.(presynPops(ii)) = presynPartners;
        else
            presynPartners = arrayfun(@(x) numel(find(model.model.fixed_variables.([char(cellToInspect) '_' char(presynPops(ii)) '_' currentType '_netcon'])(:,x))),...
                1:size(model.([cellToInspect '_V']),2),'uniformoutput',false);
            PresynapticPartnerCount.(presynPops(ii)) = presynPartners;
        end
    end
    
    % get total synaptic input
    % ... but first get time windows of interest
    timeWindows = [model.time(1) cutStartTime]; % [start stop] in ms
    tempAnalysis = analyzeSparseSim(model,'test');
    popIdx = strcmp(tempAnalysis.Population, cellToInspect);
    timeWindows = [timeWindows; cell2mat(cellfun(@(x) [x.PeakTime x.EndTime], tempAnalysis.VPMEventData{popIdx,1}, 'uniformoutput',0)')];
    
    synapticInputs.faulty_not_spiking = []; synapticInputs.faulty_spiking = []; synapticInputs.healthy = [];
    mem_pot.faulty_not_spiking = []; mem_pot.faulty_spiking = []; mem_pot.healthy = [];
    first_spikes.faulty_not_spiking = []; first_spikes.faulty_spiking = []; first_spikes.healthy = [];
    STA.faulty_not_spiking = []; STA.faulty_spiking = []; STA.healthy = [];

    for ii = 1:size(model.([cellToInspect '_V']),2)
        if ismember(ii,cellBehaviour.(cellToInspect).fns)
            synapticInputs.faulty_not_spiking{end+1} = getTotalSynInput(model,cellToInspect,ii,timeWindows);

            % get vm
            mem_pot.faulty_not_spiking{end+1} = arrayfun(@(x) model.([cellToInspect '_V'])([find(model.time == timeWindows(x,1)):find(model.time == timeWindows(x,2))], ii), 2:size(timeWindows, 1), 'UniformOutput',0);
            
            % 1st spike
            first_spikes.faulty_not_spiking{end+1} = arrayfun(@(x) tempAnalysis.VPMEventData{popIdx,1}{1, x}.FirstSpikeLatency{1, ii}, 1:length(tempAnalysis.VPMEventData{popIdx,1}), 'UniformOutput',0);

            % STA
            areThereSpikes = cell2mat(cellfun(@(x) ~isempty(x), first_spikes.faulty_not_spiking{end}, 'UniformOutput',0));
            STA_signals = [];
            if any(areThereSpikes)
                for spk_idx = 1:length(areThereSpikes)
                    if areThereSpikes(spk_idx)
                        STA_signals{spk_idx} = model.([cellToInspect '_V'])([find(round(model.time,2) == round(tempAnalysis.VPMEventData{popIdx,1}{1, spk_idx}.PeakTime + tempAnalysis.VPMEventData{popIdx,1}{1, spk_idx}.FirstSpikeLatency{1, ii} - 60,2)) :...
                            find(round(model.time,2) == round(tempAnalysis.VPMEventData{popIdx,1}{1, spk_idx}.PeakTime + tempAnalysis.VPMEventData{popIdx,1}{1, spk_idx}.FirstSpikeLatency{1, ii},2))], ii);
                    else
                        STA_signals{spk_idx} = [];
                    end
                end
                STA.faulty_not_spiking{end+1} = STA_signals;
            else
                STA.faulty_not_spiking{end+1} = repmat({[]}, 1, length(areThereSpikes));
            end

        elseif ismember(ii,cellBehaviour.(cellToInspect).fs)
            synapticInputs.faulty_spiking{end+1} = getTotalSynInput(model,cellToInspect,ii,timeWindows);

            % get vm
            mem_pot.faulty_spiking{end+1} = arrayfun(@(x) model.([cellToInspect '_V'])([find(model.time == timeWindows(x,1)):find(model.time == timeWindows(x,2))], ii), 1:size(timeWindows, 1), 'UniformOutput',0);

            % 1st spike
            first_spikes.faulty_spiking{end+1} = arrayfun(@(x) tempAnalysis.VPMEventData{popIdx,1}{1, x}.FirstSpikeLatency{1, ii}, 1:length(tempAnalysis.VPMEventData{popIdx,1}), 'UniformOutput',0);

            % STA
            areThereSpikes = cell2mat(cellfun(@(x) ~isempty(x), first_spikes.faulty_spiking{end}, 'UniformOutput',0));
            STA_signals = [];
            if any(areThereSpikes)
                for spk_idx = 1:length(areThereSpikes)
                    if areThereSpikes(spk_idx)
                        STA_signals{spk_idx} = model.([cellToInspect '_V'])([find(round(model.time,2) == round(tempAnalysis.VPMEventData{popIdx,1}{1, spk_idx}.PeakTime + tempAnalysis.VPMEventData{popIdx,1}{1, spk_idx}.FirstSpikeLatency{1, ii} - 60,2)) :...
                            find(round(model.time,2) == round(tempAnalysis.VPMEventData{popIdx,1}{1, spk_idx}.PeakTime + tempAnalysis.VPMEventData{popIdx,1}{1, spk_idx}.FirstSpikeLatency{1, ii},2))], ii);
                    else
                        STA_signals{spk_idx} = [];
                    end
                end
                STA.faulty_spiking{end+1} = STA_signals;
            else
                STA.faulty_spiking{end+1} = repmat({[]}, 1, length(areThereSpikes));
            end

        elseif ismember(ii,cellBehaviour.(cellToInspect).ht)
            synapticInputs.healthy{end+1} = getTotalSynInput(model,cellToInspect,ii,timeWindows);

            % get vm
            mem_pot.healthy{end+1} = arrayfun(@(x) model.([cellToInspect '_V'])([find(model.time == timeWindows(x,1)):find(model.time == timeWindows(x,2))], ii), 1:size(timeWindows, 1), 'UniformOutput',0);

            % 1st spike
            first_spikes.healthy{end+1} = arrayfun(@(x) tempAnalysis.VPMEventData{popIdx,1}{1, x}.FirstSpikeLatency{1, ii}, 1:length(tempAnalysis.VPMEventData{popIdx,1}), 'UniformOutput',0);

            % STA
            areThereSpikes = cell2mat(cellfun(@(x) ~isempty(x), first_spikes.healthy{end}, 'UniformOutput',0));
            STA_signals = [];
            if any(areThereSpikes)
                for spk_idx = 1:length(areThereSpikes)
                    if areThereSpikes(spk_idx)
                        STA_signals{spk_idx} = model.([cellToInspect '_V'])([find(round(model.time,2) == round(tempAnalysis.VPMEventData{popIdx,1}{1, spk_idx}.PeakTime + tempAnalysis.VPMEventData{popIdx,1}{1, spk_idx}.FirstSpikeLatency{1, ii} - 60,2)) :...
                            find(round(model.time,2) == round(tempAnalysis.VPMEventData{popIdx,1}{1, spk_idx}.PeakTime + tempAnalysis.VPMEventData{popIdx,1}{1, spk_idx}.FirstSpikeLatency{1, ii},2))], ii);
                    else
                        STA_signals{spk_idx} = [];
                    end
                end
                STA.healthy{end+1} = STA_signals;
            else
                STA.healthy{end+1} = repmat({[]}, 1, length(areThereSpikes));
            end
        end
    end
    
    % sort counted connections in categories: 
    presyn_ExcPrevLayer =[]; presyn_ExcSameLayer=[]; presyn_InhSameLayer=[]; presyn_InhPrevLayer=[];

    presyn_ExcSameLayer.faulty_not_spiking = cell2mat(PresynapticPartnerCount.(presynPops(1))(faulty_not_spiking));
    presyn_ExcSameLayer.faulty_spiking = cell2mat(PresynapticPartnerCount.(presynPops(1))(faulty_spiking));
    presyn_ExcSameLayer.healthy = cell2mat(PresynapticPartnerCount.(presynPops(1))(healthy));

    presyn_InhSameLayer.faulty_not_spiking = cell2mat(PresynapticPartnerCount.(presynPops(2))(faulty_not_spiking));
    presyn_InhSameLayer.faulty_spiking = cell2mat(PresynapticPartnerCount.(presynPops(2))(faulty_spiking));
    presyn_InhSameLayer.healthy = cell2mat(PresynapticPartnerCount.(presynPops(2))(healthy));

    presyn_ExcPrevLayer.faulty_not_spiking = cell2mat(PresynapticPartnerCount.(presynPops(3))(faulty_not_spiking));
    presyn_ExcPrevLayer.faulty_spiking = cell2mat(PresynapticPartnerCount.(presynPops(3))(faulty_spiking));
    presyn_ExcPrevLayer.healthy = cell2mat(PresynapticPartnerCount.(presynPops(3))(healthy));  
    
    % for plotting purposes
    faulty_ns_means = [mean(presyn_ExcPrevLayer.faulty_not_spiking) mean(presyn_ExcSameLayer.faulty_not_spiking) mean(presyn_InhSameLayer.faulty_not_spiking)];
    faulty_ns_stds = [std(presyn_ExcPrevLayer.faulty_not_spiking) std(presyn_ExcSameLayer.faulty_not_spiking) std(presyn_InhSameLayer.faulty_not_spiking)];

    faulty_s_means = [mean(presyn_ExcPrevLayer.faulty_spiking) mean(presyn_ExcSameLayer.faulty_spiking) mean(presyn_InhSameLayer.faulty_spiking)];
    faulty_s_stds = [std(presyn_ExcPrevLayer.faulty_spiking) std(presyn_ExcSameLayer.faulty_spiking) std(presyn_InhSameLayer.faulty_spiking)];
    
    healthy_stds = [std(presyn_ExcPrevLayer.healthy) std(presyn_ExcSameLayer.healthy) std(presyn_InhSameLayer.healthy)];
    healthy_means = [mean(presyn_ExcPrevLayer.healthy) mean(presyn_ExcSameLayer.healthy) mean(presyn_InhSameLayer.healthy)];
    
    % if there is 4 presynaptic populations, extend means and stds
    if length(presynPops)==4
        presyn_InhPrevLayer.faulty_not_spiking = cell2mat(PresynapticPartnerCount.(presynPops(4))(faulty_not_spiking));
        presyn_InhPrevLayer.faulty_spiking = cell2mat(PresynapticPartnerCount.(presynPops(4))(faulty_spiking));
        presyn_InhPrevLayer.healthy = cell2mat(PresynapticPartnerCount.(presynPops(4))(healthy));
        
        faulty_ns_means = [faulty_ns_means mean(presyn_InhPrevLayer.faulty_not_spiking)];
        faulty_ns_stds = [faulty_ns_stds mean(presyn_InhPrevLayer.faulty_not_spiking)];
        
        faulty_s_means = [faulty_s_means mean(presyn_InhPrevLayer.faulty_spiking)];
        faulty_s_stds = [faulty_s_stds mean(presyn_InhPrevLayer.faulty_spiking)];
        
        healthy_stds = [healthy_stds mean(presyn_InhPrevLayer.healthy)];
        healthy_means = [healthy_means mean(presyn_InhPrevLayer.healthy)];
    else
        presyn_InhPrevLayer.faulty_not_spiking = [];
        presyn_InhPrevLayer.faulty_spiking = [];
        presyn_InhPrevLayer.healthy = [];
    end
    
    % get an example cell per type
    if ~isempty(cellBehaviour.(cellToInspect).fns)
        example_cell=randi(length(cellBehaviour.(cellToInspect).fns),1,1);
        ec_fns = model.([cellToInspect '_V'])(:,cellBehaviour.(cellToInspect).fns(example_cell));
    else
        ec_fns = [];
    end
    if ~isempty(cellBehaviour.(cellToInspect).fs)
        example_cell=randi(length(cellBehaviour.(cellToInspect).fs),1,1);
        ec_fs = model.([cellToInspect '_V'])(:,cellBehaviour.(cellToInspect).fs(example_cell));
    else
        ec_fs = [];
    end
    if ~isempty(cellBehaviour.(cellToInspect).ht)
        example_cell=randi(length(cellBehaviour.(cellToInspect).ht),1,1);
        ec_ht = model.([cellToInspect '_V'])(:,cellBehaviour.(cellToInspect).ht(example_cell));
    else
        ec_ht = [];
    end
    
    % get cell number per type
    fns_num = length(cellBehaviour.(cellToInspect).fns);
    fs_num = length(cellBehaviour.(cellToInspect).fs);
    ht_num = length(cellBehaviour.(cellToInspect).ht);
    
    % get number of cells without presynaptic inhibition
    faulty_ns_0presynInh1 = numel(find(presyn_InhSameLayer.faulty_not_spiking==0));
    faulty_s_0presynInh1 = numel(find(presyn_InhSameLayer.faulty_spiking==0));
    healthy_0presynInh1 = numel(find(presyn_InhSameLayer.healthy==0));
    
    if length(presynPops)==4
        faulty_ns_0presynInh2 = numel(find(presyn_InhPrevLayer.faulty_not_spiking==0));
        faulty_s_0presynInh2 = numel(find(presyn_InhPrevLayer.faulty_spiking==0));
        healthy_0presynInh2 = numel(find(presyn_InhPrevLayer.healthy==0));
    else
        faulty_ns_0presynInh2=[];
        faulty_s_0presynInh2=[];
        healthy_0presynInh2=[];
    end
    
    % total syn input sorted
    if ~isempty(synapticInputs.healthy)
        tsyn_ht_events_ExcPrevLayer = cell2mat(... % events
            cellfun(@(x) IntegrateEventsSeparately_PerNeuron(x, presynPops(3)), synapticInputs.healthy, 'uniformoutput',false)...
            );

        tsyn_ht_full_ExcPrevLayer = cell2mat(... % simulation without settling time
            arrayfun(@(x) trapz(synapticInputs.healthy{1, x}.(presynPops(3)).fullSimulation{1, 1}),... % full simulation
            1:length(synapticInputs.healthy),'uniformoutput',false)...
            ) - cell2mat(...
            arrayfun(@(x) trapz(synapticInputs.healthy{1, x}.(presynPops(3)).timeWindows{1, 1}),... % settling time
            1:length(synapticInputs.healthy),'uniformoutput',false)...
            );

        tsyn_ht_events_ExcSameLayer = cell2mat(... % events
            cellfun(@(x) IntegrateEventsSeparately_PerNeuron(x, presynPops(1)), synapticInputs.healthy, 'uniformoutput',false)...
            );

        tsyn_ht_full_ExcSameLayer = cell2mat(... % simulation without settling time
            arrayfun(@(x) trapz(synapticInputs.healthy{1, x}.(presynPops(1)).fullSimulation{1, 1}),... % full simulation
            1:length(synapticInputs.healthy),'uniformoutput',false)...
            ) - cell2mat(...
            arrayfun(@(x) trapz(synapticInputs.healthy{1, x}.(presynPops(1)).timeWindows{1, 1}),... % settling time
            1:length(synapticInputs.healthy),'uniformoutput',false)...
            );

        tsyn_ht_events_InhSameLayer = cell2mat(... % events
            cellfun(@(x) IntegrateEventsSeparately_PerNeuron(x, presynPops(2)), synapticInputs.healthy, 'uniformoutput',false)...
            );

        tsyn_ht_full_InhSameLayer = cell2mat(... % simulation without settling time
            arrayfun(@(x) trapz(synapticInputs.healthy{1, x}.(presynPops(2)).fullSimulation{1, 1}),... % full simulation
            1:length(synapticInputs.healthy),'uniformoutput',false)...
            ) - cell2mat(...
            arrayfun(@(x) trapz(synapticInputs.healthy{1, x}.(presynPops(2)).timeWindows{1, 1}),... % settling time
            1:length(synapticInputs.healthy),'uniformoutput',false)...
            );
        
        if length(presynPops)==4
            tsyn_ht_events_InhPrevLayer = cell2mat(... % events
                cellfun(@(x) IntegrateEventsSeparately_PerNeuron(x, presynPops(4)), synapticInputs.healthy, 'uniformoutput',false)...
                );

            tsyn_ht_full_InhPrevLayer = cell2mat(... % simulation without settling time
                arrayfun(@(x) trapz(synapticInputs.healthy{1, x}.(presynPops(4)).fullSimulation{1, 1}),... % full simulation
                1:length(synapticInputs.healthy),'uniformoutput',false)...
                ) - cell2mat(...
                arrayfun(@(x) trapz(synapticInputs.healthy{1, x}.(presynPops(4)).timeWindows{1, 1}),... % settling time
                1:length(synapticInputs.healthy),'uniformoutput',false)...
                );
        else
            tsyn_ht_events_InhPrevLayer = [];
            tsyn_ht_full_InhPrevLayer = [];
        end
    else
        tsyn_ht_events_ExcPrevLayer = [];
        tsyn_ht_full_ExcPrevLayer = [];
        tsyn_ht_events_ExcSameLayer = [];
        tsyn_ht_full_ExcSameLayer = [];
        tsyn_ht_events_InhSameLayer = [];
        tsyn_ht_full_InhSameLayer = [];
        tsyn_ht_events_InhPrevLayer = [];
        tsyn_ht_full_InhPrevLayer = [];
    end
    
    if ~isempty(synapticInputs.faulty_not_spiking)
        tsyn_fns_events_ExcPrevLayer = cell2mat(... % events
            cellfun(@(x) IntegrateEventsSeparately_PerNeuron(x, presynPops(3)), synapticInputs.faulty_not_spiking, 'uniformoutput',false)...
            );

        tsyn_fns_full_ExcPrevLayer = cell2mat(... % simulation without settling time
            arrayfun(@(x) trapz(synapticInputs.faulty_not_spiking{1, x}.(presynPops(3)).fullSimulation{1, 1}),... % full simulation
            1:length(synapticInputs.faulty_not_spiking),'uniformoutput',false)...
            ) - cell2mat(...
            arrayfun(@(x) trapz(synapticInputs.faulty_not_spiking{1, x}.(presynPops(3)).timeWindows{1, 1}),... % settling time
            1:length(synapticInputs.faulty_not_spiking),'uniformoutput',false)...
            );
        
        tsyn_fns_events_ExcSameLayer = cell2mat(... % events
            cellfun(@(x) IntegrateEventsSeparately_PerNeuron(x, presynPops(1)), synapticInputs.faulty_not_spiking, 'uniformoutput',false)...
            );

        tsyn_fns_full_ExcSameLayer = cell2mat(... % simulation without settling time
            arrayfun(@(x) trapz(synapticInputs.faulty_not_spiking{1, x}.(presynPops(1)).fullSimulation{1, 1}),... % full simulation
            1:length(synapticInputs.faulty_not_spiking),'uniformoutput',false)...
            ) - cell2mat(...
            arrayfun(@(x) trapz(synapticInputs.faulty_not_spiking{1, x}.(presynPops(1)).timeWindows{1, 1}),... % settling time
            1:length(synapticInputs.faulty_not_spiking),'uniformoutput',false)...
            );
        
        tsyn_fns_events_InhSameLayer = cell2mat(... % events
            cellfun(@(x) IntegrateEventsSeparately_PerNeuron(x, presynPops(2)), synapticInputs.faulty_not_spiking, 'uniformoutput',false)...
            );

        tsyn_fns_full_InhSameLayer = cell2mat(... % simulation without settling time
            arrayfun(@(x) trapz(synapticInputs.faulty_not_spiking{1, x}.(presynPops(2)).fullSimulation{1, 1}),... % full simulation
            1:length(synapticInputs.faulty_not_spiking),'uniformoutput',false)...
            ) - cell2mat(...
            arrayfun(@(x) trapz(synapticInputs.faulty_not_spiking{1, x}.(presynPops(2)).timeWindows{1, 1}),... % settling time
            1:length(synapticInputs.faulty_not_spiking),'uniformoutput',false)...
            );
        
        if length(presynPops)==4
            tsyn_fns_events_InhPrevLayer = cell2mat(... % events
                cellfun(@(x) IntegrateEventsSeparately_PerNeuron(x, presynPops(4)), synapticInputs.faulty_not_spiking, 'uniformoutput',false)...
                );

            tsyn_fns_full_InhPrevLayer = cell2mat(... % simulation without settling time
                arrayfun(@(x) trapz(synapticInputs.faulty_not_spiking{1, x}.(presynPops(4)).fullSimulation{1, 1}),... % full simulation
                1:length(synapticInputs.faulty_not_spiking),'uniformoutput',false)...
                ) - cell2mat(...
                arrayfun(@(x) trapz(synapticInputs.faulty_not_spiking{1, x}.(presynPops(4)).timeWindows{1, 1}),... % settling time
                1:length(synapticInputs.faulty_not_spiking),'uniformoutput',false)...
                );
        else
            tsyn_fns_events_InhPrevLayer = [];
            tsyn_fns_full_InhPrevLayer = [];
        end
    else
        tsyn_fns_events_ExcPrevLayer = [];
        tsyn_fns_full_ExcPrevLayer = [];
        tsyn_fns_events_ExcSameLayer = [];
        tsyn_fns_full_ExcSameLayer = [];
        tsyn_fns_events_InhSameLayer = [];
        tsyn_fns_full_InhSameLayer = [];
        tsyn_fns_events_InhPrevLayer = [];
        tsyn_fns_full_InhPrevLayer = [];
    end
    
    if ~isempty(synapticInputs.faulty_spiking)
        tsyn_fs_events_ExcPrevLayer = cell2mat(... % events
            cellfun(@(x) IntegrateEventsSeparately_PerNeuron(x, presynPops(3)), synapticInputs.faulty_spiking, 'uniformoutput',false)...
            );

        tsyn_fs_full_ExcPrevLayer = cell2mat(... % simulation without settling time
            arrayfun(@(x) trapz(synapticInputs.faulty_spiking{1, x}.(presynPops(3)).fullSimulation{1, 1}),... % full simulation
            1:length(synapticInputs.faulty_spiking),'uniformoutput',false)...
            ) - cell2mat(...
            arrayfun(@(x) trapz(synapticInputs.faulty_spiking{1, x}.(presynPops(3)).timeWindows{1, 1}),... % settling time
            1:length(synapticInputs.faulty_spiking),'uniformoutput',false)...
            );
        
        tsyn_fs_events_ExcSameLayer = cell2mat(... % events
            cellfun(@(x) IntegrateEventsSeparately_PerNeuron(x, presynPops(1)), synapticInputs.faulty_spiking, 'uniformoutput',false)...
            );

        tsyn_fs_full_ExcSameLayer = cell2mat(... % simulation without settling time
            arrayfun(@(x) trapz(synapticInputs.faulty_spiking{1, x}.(presynPops(1)).fullSimulation{1, 1}),... % full simulation
            1:length(synapticInputs.faulty_spiking),'uniformoutput',false)...
            ) - cell2mat(...
            arrayfun(@(x) trapz(synapticInputs.faulty_spiking{1, x}.(presynPops(1)).timeWindows{1, 1}),... % settling time
            1:length(synapticInputs.faulty_spiking),'uniformoutput',false)...
            );
        
        tsyn_fs_events_InhSameLayer = cell2mat(... % events
            cellfun(@(x) IntegrateEventsSeparately_PerNeuron(x, presynPops(2)), synapticInputs.faulty_spiking, 'uniformoutput',false)...
            );

        tsyn_fs_full_InhSameLayer = cell2mat(... % simulation without settling time
            arrayfun(@(x) trapz(synapticInputs.faulty_spiking{1, x}.(presynPops(2)).fullSimulation{1, 1}),... % full simulation
            1:length(synapticInputs.faulty_spiking),'uniformoutput',false)...
            ) - cell2mat(...
            arrayfun(@(x) trapz(synapticInputs.faulty_spiking{1, x}.(presynPops(2)).timeWindows{1, 1}),... % settling time
            1:length(synapticInputs.faulty_spiking),'uniformoutput',false)...
            );
        
        if length(presynPops)==4
            tsyn_fs_events_InhPrevLayer = cell2mat(... % events
                cellfun(@(x) IntegrateEventsSeparately_PerNeuron(x, presynPops(4)), synapticInputs.faulty_spiking, 'uniformoutput',false)...
                );

            tsyn_fs_full_InhPrevLayer = cell2mat(... % simulation without settling time
                arrayfun(@(x) trapz(synapticInputs.faulty_spiking{1, x}.(presynPops(4)).fullSimulation{1, 1}),... % full simulation
                1:length(synapticInputs.faulty_spiking),'uniformoutput',false)...
                ) - cell2mat(...
                arrayfun(@(x) trapz(synapticInputs.faulty_spiking{1, x}.(presynPops(4)).timeWindows{1, 1}),... % settling time
                1:length(synapticInputs.faulty_spiking),'uniformoutput',false)...
                );
        else
            tsyn_fs_events_InhPrevLayer = [];
            tsyn_fs_full_InhPrevLayer = [];
        end
    else
        tsyn_fs_events_ExcPrevLayer = [];
        tsyn_fs_full_ExcPrevLayer = [];
        tsyn_fs_events_ExcSameLayer = [];
        tsyn_fs_full_ExcSameLayer = [];
        tsyn_fs_events_InhSameLayer = [];
        tsyn_fs_full_InhSameLayer = [];
        tsyn_fs_events_InhPrevLayer = [];
        tsyn_fs_full_InhPrevLayer = [];
    end
    
    
    % plot another figure if the user wants
    if plon(2)==1
        f1001=figure; tiledlayout(3,3); 
        nexttile; hold on;
        plot([cutStartTime cutStartTime], [-100 50], ':k')
        if ~isempty(ec_fns)
            plot(time(1):time(3):time(2),ec_fns,'k','linewidth',1.2)
        end
        title('Faulty: not spiking (FNS)')
        xlabel('time (ms)')
        set(gca,'fontsize',14)
        xlim([0 time(2)])

        nexttile; hold on;
        plot([cutStartTime cutStartTime], [-100 50], ':k')
        if ~isempty(ec_fs)
            plot(time(1):time(3):time(2),ec_fs,'k','linewidth',1.2)
        end
        title('Faulty: spiking (FS)')
        xlabel('time (ms)')
        set(gca,'fontsize',14)
        xlim([0 time(2)])

        nexttile; hold on;
        plot([cutStartTime cutStartTime], [-100 50], ':k')
        if ~isempty(ec_ht)
            plot(time(1):time(3):time(2),ec_ht,'k','linewidth',1.2)
        end
        title('Healthy (HT)')
        xlabel('time (ms)')
        set(gca,'fontsize',14)
        xlim([0 time(2)])
        
        nexttile([2 1]); hold on;
        br1 = bar([1:length(presynPops)],faulty_ns_means,'grouped','k','BarWidth',0.5,'EdgeColor','none','FaceAlpha',0.2);
        br2 = bar([1:length(presynPops)]+4,faulty_s_means,'grouped','b','BarWidth',0.5,'EdgeColor','none','FaceAlpha',0.2);
        br3 = bar([1:length(presynPops)]+8,healthy_means,'grouped','facecolor',[0.008, 0.424, 0.271],'BarWidth',0.5,'EdgeColor','none','FaceAlpha',0.3);
        errorbar([br1.XData; br2.XData; br3.XData],[faulty_ns_means;faulty_s_means;healthy_means],[faulty_ns_stds;faulty_s_stds;healthy_stds],'linestyle','none','color','k','LineWidth',1.4);
        try
            legend(['FNS=' num2str(fns_num) 'cells'],['FS=' num2str(fs_num) 'cells'],['HT=' num2str(ht_num) 'cells'],'location','bestoutside')
        catch ME
            ax = get(gca,'Children');
            bx1=findobj(ax(1).Children,'Tag','Box'); bx1=bx1(1);
            bx2=findobj(ax(2).Children,'Tag','Box'); bx2=bx2(1);
            bx3=findobj(ax(3).Children,'Tag','Box'); bx3=bx3(1);
            legend([bx1,bx2,bx3],['FNS=' num2str(fns_num) 'cells'],['FS=' num2str(fs_num) 'cells'],['HT=' num2str(ht_num) 'cells'],'location','bestoutside')
        end
        xlabel('Presynaptic neuron')
        ylabel('# connections')
        ylim([0 25])
        ax=gca;
        ax.XTickMode='manual';
        ax.XTick=[br1.XData br2.XData br3.XData];
        if length(br1.XData)==3
            set(gca,'fontsize',14,'xticklabel',[presynPops(3) presynPops(1) presynPops(2)])
        else
            set(gca,'fontsize',14,'xticklabel',[presynPops(3) presynPops(1) presynPops(2) presynPops(4)])
        end
        title({'Presynaptic connections'; ['per ' cellToInspect ' cell type']})
        
        nexttile; hold on
        br1 = bar(1,faulty_ns_0presynInh1/length(presyn_InhSameLayer.faulty_not_spiking),'k','BarWidth',0.5,'EdgeColor','none','FaceAlpha',0.2);
        br2 = bar(2,faulty_s_0presynInh1/length(presyn_InhSameLayer.faulty_spiking),'b','BarWidth',0.5,'EdgeColor','none','FaceAlpha',0.2);
        br3 = bar(3,healthy_0presynInh1/length(presyn_InhSameLayer.healthy),'facecolor',[0.008, 0.424, 0.271],'BarWidth',0.5,'EdgeColor','none','FaceAlpha',0.3);
        xlabel('Presynapse type')
        ylabel('% of cells/type')
        title({'Fraction of cells'; 'with no inh. presynapse'})    
        legend([num2str(faulty_ns_0presynInh1) '/' num2str(length(presyn_InhSameLayer.faulty_not_spiking))],...
            [num2str(faulty_s_0presynInh1) '/' num2str(length(presyn_InhSameLayer.faulty_spiking))],...
            [num2str(healthy_0presynInh1) '/' num2str(length(presyn_InhSameLayer.healthy))])
        ylim([0 1])
        ax=gca;
        ax.XTickMode='manual';
        ax.XTick=[br1.XData br2.XData br3.XData];
        set(gca,'fontsize',14,'xticklabel',presynPops(2))
        
        if length(presynPops)==4
            br4 = bar(4,faulty_ns_0presynInh2/length(presyn_InhPrevLayer.faulty_not_spiking),'k','BarWidth',0.5,'EdgeColor','none','FaceAlpha',0.4);
            br5 = bar(5,faulty_s_0presynInh2/length(presyn_InhPrevLayer.faulty_spiking),'b','BarWidth',0.5,'EdgeColor','none','FaceAlpha',0.5);
            br6 = bar(6,healthy_0presynInh2/length(presyn_InhPrevLayer.healthy),'facecolor',[0.008, 0.424, 0.271],'BarWidth',0.5,'EdgeColor','none','FaceAlpha',0.6);
            legend([num2str(faulty_ns_0presynInh1) '/' num2str(length(presyn_InhSameLayer.faulty_not_spiking))],...
                [num2str(faulty_s_0presynInh1) '/' num2str(length(presyn_InhSameLayer.faulty_spiking))],...
                [num2str(healthy_0presynInh1) '/' num2str(length(presyn_InhSameLayer.healthy))],...
                [num2str(faulty_ns_0presynInh2) '/' num2str(length(presyn_InhPrevLayer.faulty_not_spiking))],...
                [num2str(faulty_s_0presynInh2) '/' num2str(length(presyn_InhPrevLayer.faulty_spiking))],...
                [num2str(healthy_0presynInh2) '/' num2str(length(presyn_InhPrevLayer.healthy))])
            ylim([0 1])
            ax=gca;
            ax.XTickMode='manual';
            ax.XTick=[br1.XData br2.XData br3.XData br4.XData br5.XData br6.XData];
            set(gca,'fontsize',14,'xticklabel',[repelem(presynPops(2),3) repelem(presynPops(4),3)])
        end
        
        nexttile; hold on; % number of excitatory partners in cells without inhibition
        fns_0presynInh1_presynExcPrevLayer_mean = mean(presyn_ExcPrevLayer.faulty_not_spiking(presyn_InhSameLayer.faulty_not_spiking==0));
        fns_0presynInh1_presynExcSameLayer_mean = mean(presyn_ExcSameLayer.faulty_not_spiking(presyn_InhSameLayer.faulty_not_spiking==0));
        fns_0presynInh1_presynExcPrevLayer_std = std(presyn_ExcPrevLayer.faulty_not_spiking(presyn_InhSameLayer.faulty_not_spiking==0));
        fns_0presynInh1_presynExcSameLayer_std = std(presyn_ExcSameLayer.faulty_not_spiking(presyn_InhSameLayer.faulty_not_spiking==0));
        ht_0presynInh1_presynExcPrevLayer_mean = mean(presyn_ExcPrevLayer.healthy(presyn_InhSameLayer.healthy==0));
        ht_0presynInh1_presynExcSameLayer_mean = mean(presyn_ExcSameLayer.healthy(presyn_InhSameLayer.healthy==0));
        ht_0presynInh1_presynExcPrevLayer_std = std(presyn_ExcPrevLayer.healthy(presyn_InhSameLayer.healthy==0));
        ht_0presynInh1_presynExcSameLayer_std = std(presyn_ExcSameLayer.healthy(presyn_InhSameLayer.healthy==0));

        fns_1presynInh1_presynExcPrevLayer_mean = mean(presyn_ExcPrevLayer.faulty_not_spiking(presyn_InhSameLayer.faulty_not_spiking~=0));
        fns_1presynInh1_presynExcSameLayer_mean = mean(presyn_ExcSameLayer.faulty_not_spiking(presyn_InhSameLayer.faulty_not_spiking~=0));
        fns_1presynInh1_presynExcPrevLayer_std = std(presyn_ExcPrevLayer.faulty_not_spiking(presyn_InhSameLayer.faulty_not_spiking~=0));
        fns_1presynInh1_presynExcSameLayer_std = std(presyn_ExcSameLayer.faulty_not_spiking(presyn_InhSameLayer.faulty_not_spiking~=0));
        ht_1presynInh1_presynExcPrevLayer_mean = mean(presyn_ExcPrevLayer.healthy(presyn_InhSameLayer.healthy~=0));
        ht_1presynInh1_presynExcSameLayer_mean = mean(presyn_ExcSameLayer.healthy(presyn_InhSameLayer.healthy~=0));
        ht_1presynInh1_presynExcPrevLayer_std = std(presyn_ExcPrevLayer.healthy(presyn_InhSameLayer.healthy~=0));
        ht_1presynInh1_presynExcSameLayer_std = std(presyn_ExcSameLayer.healthy(presyn_InhSameLayer.healthy~=0));

        br1 = bar([1 2],[fns_0presynInh1_presynExcPrevLayer_mean fns_0presynInh1_presynExcSameLayer_mean],'k','BarWidth',0.5,'EdgeColor','none','FaceAlpha',0.2);
        br2 = bar([3 4],[fns_1presynInh1_presynExcPrevLayer_mean fns_1presynInh1_presynExcSameLayer_mean],'EdgeColor','k','BarWidth',0.5,'faceColor','none','EdgeAlpha',0.8);
        br3 = bar([6 7],[ht_0presynInh1_presynExcPrevLayer_mean ht_0presynInh1_presynExcSameLayer_mean],'facecolor',[0.008, 0.424, 0.271],'BarWidth',0.5,'EdgeColor','none','FaceAlpha',0.3);
        br4 = bar([8 9],[ht_1presynInh1_presynExcPrevLayer_mean ht_1presynInh1_presynExcSameLayer_mean],'EdgeColor',[0.008, 0.424, 0.271],'BarWidth',0.5,'faceColor','none','EdgeAlpha',0.7);
        errorbar([br1.XData br2.XData; br3.XData br4.XData], [fns_0presynInh1_presynExcPrevLayer_mean fns_0presynInh1_presynExcSameLayer_mean fns_1presynInh1_presynExcPrevLayer_mean fns_1presynInh1_presynExcSameLayer_mean; ht_0presynInh1_presynExcPrevLayer_mean ht_0presynInh1_presynExcSameLayer_mean ht_1presynInh1_presynExcPrevLayer_mean ht_1presynInh1_presynExcSameLayer_mean],...
            [fns_0presynInh1_presynExcPrevLayer_std fns_0presynInh1_presynExcSameLayer_std fns_1presynInh1_presynExcPrevLayer_std fns_1presynInh1_presynExcSameLayer_std; ht_0presynInh1_presynExcPrevLayer_std ht_0presynInh1_presynExcSameLayer_std ht_1presynInh1_presynExcPrevLayer_std ht_1presynInh1_presynExcSameLayer_std],'linestyle','none','color','k','LineWidth',1.4)
        ax=gca;
        ax.XTickMode='manual';
        ax.XTick=[br1.XData br2.XData br3.XData br4.XData];
        set(gca,'fontsize',14,'xticklabel',presynPops([3 1]))
        legend(sprintf('FNS -%s', presynPops(2)), sprintf('FNS +%s', presynPops(2)), sprintf('HT -%s', presynPops(2)), sprintf('HT +%s', presynPops(2)), 'location', 'eastoutside')
        xlabel('Presynapse type')
        ylabel('# partners')
        title({'Excitatory presynaptic partner count'; 'in cells with (+) and without (-) inh. presynapse'})
        
        if length(presynPops)==4
            fns_0presynInh2_presynExcPrevLayer_mean = mean(presyn_ExcPrevLayer.faulty_not_spiking(presyn_InhPrevLayer.faulty_not_spiking==0));
            fns_0presynInh2_presynExcSameLayer_mean = mean(presyn_ExcSameLayer.faulty_not_spiking(presyn_InhPrevLayer.faulty_not_spiking==0));
            fns_0presynInh2_presynExcPrevLayer_std = std(presyn_ExcPrevLayer.faulty_not_spiking(presyn_InhPrevLayer.faulty_not_spiking==0));
            fns_0presynInh2_presynExcSameLayer_std = std(presyn_ExcSameLayer.faulty_not_spiking(presyn_InhPrevLayer.faulty_not_spiking==0));
            ht_0presynInh2_presynExcPrevLayer_mean = mean(presyn_ExcPrevLayer.healthy(presyn_InhPrevLayer.healthy==0));
            ht_0presynInh2_presynExcSameLayer_mean = mean(presyn_ExcSameLayer.healthy(presyn_InhPrevLayer.healthy==0));
            ht_0presynInh2_presynExcPrevLayer_std = std(presyn_ExcPrevLayer.healthy(presyn_InhPrevLayer.healthy==0));
            ht_0presynInh2_presynExcSameLayer_std = std(presyn_ExcSameLayer.healthy(presyn_InhPrevLayer.healthy==0));

            fns_1presynInh2_presynExcPrevLayer_mean = mean(presyn_ExcPrevLayer.faulty_not_spiking(presyn_InhPrevLayer.faulty_not_spiking~=0));
            fns_1presynInh2_presynExcSameLayer_mean = mean(presyn_ExcSameLayer.faulty_not_spiking(presyn_InhPrevLayer.faulty_not_spiking~=0));
            fns_1presynInh2_presynExcPrevLayer_std = std(presyn_ExcPrevLayer.faulty_not_spiking(presyn_InhPrevLayer.faulty_not_spiking~=0));
            fns_1presynInh2_presynExcSameLayer_std = std(presyn_ExcSameLayer.faulty_not_spiking(presyn_InhPrevLayer.faulty_not_spiking~=0));
            ht_1presynInh2_presynExcPrevLayer_mean = mean(presyn_ExcPrevLayer.healthy(presyn_InhPrevLayer.healthy~=0));
            ht_1presynInh2_presynExcSameLayer_mean = mean(presyn_ExcSameLayer.healthy(presyn_InhPrevLayer.healthy~=0));
            ht_1presynInh2_presynExcPrevLayer_std = std(presyn_ExcPrevLayer.healthy(presyn_InhPrevLayer.healthy~=0));
            ht_1presynInh2_presynExcSameLayer_std = std(presyn_ExcSameLayer.healthy(presyn_InhPrevLayer.healthy~=0));
       
            br5 = bar([1 2]+10,[fns_0presynInh2_presynExcPrevLayer_mean fns_0presynInh2_presynExcSameLayer_mean],'k','BarWidth',0.5,'EdgeColor','none','FaceAlpha',0.4);
            br6 = bar([3 4]+10,[fns_1presynInh2_presynExcPrevLayer_mean fns_1presynInh2_presynExcSameLayer_mean],'EdgeColor','k','BarWidth',0.5,'faceColor','none','EdgeAlpha',1);
            br7 = bar([6 7]+10,[ht_0presynInh2_presynExcPrevLayer_mean ht_0presynInh2_presynExcSameLayer_mean],'facecolor',[0.008, 0.424, 0.271],'BarWidth',0.5,'EdgeColor','none','FaceAlpha',0.6);
            br8 = bar([8 9]+10,[ht_1presynInh2_presynExcPrevLayer_mean ht_1presynInh2_presynExcSameLayer_mean],'EdgeColor',[0.008, 0.424, 0.271],'BarWidth',0.5,'faceColor','none','EdgeAlpha',1);
            errorbar([br5.XData br6.XData; br7.XData br8.XData], [fns_0presynInh2_presynExcPrevLayer_mean fns_0presynInh2_presynExcSameLayer_mean fns_1presynInh2_presynExcPrevLayer_mean fns_1presynInh2_presynExcSameLayer_mean; ht_0presynInh2_presynExcPrevLayer_mean ht_0presynInh2_presynExcSameLayer_mean ht_1presynInh2_presynExcPrevLayer_mean ht_1presynInh2_presynExcSameLayer_mean],...
                [fns_0presynInh2_presynExcPrevLayer_std fns_0presynInh2_presynExcSameLayer_std fns_1presynInh2_presynExcPrevLayer_std fns_1presynInh2_presynExcSameLayer_std; ht_0presynInh2_presynExcPrevLayer_std ht_0presynInh2_presynExcSameLayer_std ht_1presynInh2_presynExcPrevLayer_std ht_1presynInh2_presynExcSameLayer_std],'linestyle','none','color','k','LineWidth',1.2)
            
            ax.XTick=[br1.XData br2.XData br3.XData br4.XData br5.XData br6.XData br7.XData br8.XData];
            lgd=get(gca,'Children');
            legend(flipud([lgd(5:end-8); lgd(end-3:end)]),sprintf('FNS -%s', presynPops(2)), sprintf('FNS +%s', presynPops(2)), sprintf('HT -%s', presynPops(2)), sprintf('HT +%s', presynPops(2)),...
                sprintf('FNS -%s', presynPops(4)), sprintf('FNS +%s', presynPops(4)), sprintf('HT -%s', presynPops(4)), sprintf('HT +%s', presynPops(4)),'location', 'eastoutside')   
        end
        
        nexttile; hold on;
        if length(presynPops)==4
            % settling time
            try % if there is only 1 fns, then boxplot gets confused and it groups all three tsyn as one, so then "'positions', [1 2 3]" throws an error
                boxplot([tsyn_fns_events_ExcPrevLayer; tsyn_fns_events_ExcSameLayer; tsyn_fns_events_InhSameLayer; tsyn_fns_events_InhPrevLayer]','positions',[1 2 3 4],'plotstyle','compact','boxstyle','filled','colors','k','medianstyle','line')
            catch ME
                boxplot([[tsyn_fns_events_ExcPrevLayer; tsyn_fns_events_ExcSameLayer; tsyn_fns_events_InhSameLayer; tsyn_fns_events_InhPrevLayer]'; nan(1,4)],'positions',[1 2 3 4],'plotstyle','compact','boxstyle','filled','colors','k','medianstyle','line')
            end
            try
                boxplot([tsyn_fs_events_ExcPrevLayer; tsyn_fs_events_ExcSameLayer; tsyn_fs_events_InhSameLayer; tsyn_fs_events_InhPrevLayer]','positions',[1 2 3 4]+5,'plotstyle','compact','boxstyle','filled','colors','k','medianstyle','line')
            catch ME
                boxplot([[tsyn_fs_events_ExcPrevLayer; tsyn_fs_events_ExcSameLayer; tsyn_fs_events_InhSameLayer; tsyn_fs_events_InhPrevLayer]'; nan(1,4)],'positions',[1 2 3 4]+5,'plotstyle','compact','boxstyle','filled','colors','k','medianstyle','line')
            end
            boxplot([tsyn_ht_events_ExcPrevLayer; tsyn_ht_events_ExcSameLayer; tsyn_ht_events_InhSameLayer; tsyn_ht_events_InhPrevLayer]','positions',[1 2 3 4]+10,'plotstyle','compact','boxstyle','filled','colors','k','medianstyle','line')
        elseif length(presynPops)==3
            % settling time
            try % if there is only 1 fns, then boxplot gets confused and it groups all three tsyn as one, so then "'positions', [1 2 3]" throws an error
                boxplot([tsyn_fns_events_ExcPrevLayer; tsyn_fns_events_ExcSameLayer; tsyn_fns_events_InhSameLayer]','positions',[1 2 3],'plotstyle','compact','boxstyle','filled','colors','k','medianstyle','line')
            catch ME
                boxplot([[tsyn_fns_events_ExcPrevLayer; tsyn_fns_events_ExcSameLayer; tsyn_fns_events_InhSameLayer]'; nan(1,3)],'positions',[1 2 3],'plotstyle','compact','boxstyle','filled','colors','k','medianstyle','line')
            end
            try
                boxplot([tsyn_fs_events_ExcPrevLayer; tsyn_fs_events_ExcSameLayer; tsyn_fs_events_InhSameLayer]','positions',[1 2 3]+4,'plotstyle','compact','boxstyle','filled','colors','k','medianstyle','line')
            catch ME
                boxplot([[tsyn_fs_events_ExcPrevLayer; tsyn_fs_events_ExcSameLayer; tsyn_fs_events_InhSameLayer]'; nan(1,3)],'positions',[1 2 3]+4,'plotstyle','compact','boxstyle','filled','colors','k','medianstyle','line')
            end
            boxplot([tsyn_ht_events_ExcPrevLayer; tsyn_ht_events_ExcSameLayer; tsyn_ht_events_InhSameLayer]','positions',[1 2 3]+8,'plotstyle','compact','boxstyle','filled','colors','k','medianstyle','line')
        end
        
        xlabel('Presynapse type')
        ylabel('Current (pA)')
        ylims = abs(ylim);
        ylim([-max(ylims) max(ylims)])
        plot(xlim, [0 0], ':k')
        ax=gca;
        ax.XTickMode='manual';
        if length(presynPops)==3
            ax.XTick=[1:11];
            set(gca,'fontsize',14,'xticklabel',[presynPops(3) presynPops(1) presynPops(2) " "])
            h = findobj(gca,'Tag','Box');
            for jj=1:length(h)
                if jj<4
                    h(jj).Color = [0.008, 0.424, 0.271 0.3]; % ht
                elseif jj>3 && jj<7
                    h(jj).Color = [0 0 1 0.3]; % fs
                elseif jj>6
                    h(jj).Color = [0 0 0 0.3]; % fns
                end
            end
        else
            ax.XTick=[1:14];
            set(gca,'fontsize',14,'xticklabel',[presynPops(3) presynPops(1) presynPops(2) presynPops(4) " "])
            h = findobj(gca,'Tag','Box');
            for jj=1:length(h)
                if jj<5
                    h(jj).Color = [0.008, 0.424, 0.271 0.3]; % ht
                elseif jj>4 && jj<8
                    h(jj).Color = [0 0 1 0.3]; % fs
                elseif jj>7
                    h(jj).Color = [0 0 0 0.3]; % fns
                end
            end
        end
        title({'Total synaptic input (Events)'; ['per ' char(cellToInspect) ' cell']})

        % simulation time
        nexttile; hold on;
        if length(presynPops)==3
            try
                boxplot([tsyn_fns_full_ExcPrevLayer; tsyn_fns_full_ExcSameLayer; tsyn_fns_full_InhSameLayer]','positions',[1 2 3],'plotstyle','compact','boxstyle','filled','colors','k', 'medianstyle','line')
            catch ME
                boxplot([[tsyn_fns_full_ExcPrevLayer; tsyn_fns_full_ExcSameLayer; tsyn_fns_full_InhSameLayer]'; nan(1,3)],'positions',[1 2 3],'plotstyle','compact','boxstyle','filled','colors','k', 'medianstyle','line')
            end
            try
                boxplot([tsyn_fs_full_ExcPrevLayer; tsyn_fs_full_ExcSameLayer; tsyn_fs_full_InhSameLayer]','positions',[1 2 3]+4,'plotstyle','compact','boxstyle','filled','colors','k','medianstyle','line')
            catch ME
                boxplot([[tsyn_fs_full_ExcPrevLayer; tsyn_fs_full_ExcSameLayer; tsyn_fs_full_InhSameLayer]'; nan(1,3)],'positions',[1 2 3]+4,'plotstyle','compact','boxstyle','filled','colors','k','medianstyle','line')
            end
            boxplot([tsyn_ht_full_ExcPrevLayer; tsyn_ht_full_ExcSameLayer; tsyn_ht_full_InhSameLayer]','positions',[1 2 3]+8,'plotstyle','compact','boxstyle','filled','colors','k','medianstyle','line')
        else
            try
                boxplot([tsyn_fns_full_ExcPrevLayer; tsyn_fns_full_ExcSameLayer; tsyn_fns_full_InhSameLayer; tsyn_fns_full_InhPrevLayer]','positions',[1 2 3 4],'plotstyle','compact','boxstyle','filled','colors','k', 'medianstyle','line')
            catch ME
                boxplot([[tsyn_fns_full_ExcPrevLayer; tsyn_fns_full_ExcSameLayer; tsyn_fns_full_InhSameLayer; tsyn_fns_full_InhPrevLayer]'; nan(1,4)],'positions',[1 2 3 4],'plotstyle','compact','boxstyle','filled','colors','k', 'medianstyle','line')
            end
            try
                boxplot([tsyn_fs_full_ExcPrevLayer; tsyn_fs_full_ExcSameLayer; tsyn_fs_full_InhSameLayer; tsyn_fs_full_InhPrevLayer]','positions',[1 2 3 4]+5,'plotstyle','compact','boxstyle','filled','colors','k','medianstyle','line')
            catch ME
                boxplot([[tsyn_fs_full_ExcPrevLayer; tsyn_fs_full_ExcSameLayer; tsyn_fs_full_InhSameLayer; tsyn_fs_full_InhPrevLayer]'; nan(1,4)],'positions',[1 2 3 4]+5,'plotstyle','compact','boxstyle','filled','colors','k','medianstyle','line')
            end
            boxplot([tsyn_ht_full_ExcPrevLayer; tsyn_ht_full_ExcSameLayer; tsyn_ht_full_InhSameLayer; tsyn_ht_full_InhPrevLayer]','positions',[1 2 3 4]+10,'plotstyle','compact','boxstyle','filled','colors','k','medianstyle','line')
        end
        xlabel('Presynapse type')
        ylabel('Current (pA)')
        ylims = abs(ylim);
        ylim([-max(ylims) max(ylims)])
        plot(xlim, [0 0], ':k')
        ax=gca;
        ax.XTickMode='manual';
        if length(presynPops)==3
            ax.XTick=[1:11];
            set(gca,'fontsize',14,'xticklabel',[presynPops(3) presynPops(1) presynPops(2) " "])
            h = findobj(gca,'Tag','Box');
            for jj=1:length(h)
                if jj<4
                    h(jj).Color = [0.008, 0.424, 0.271 0.3]; % ht
                elseif jj>3 && jj<7
                    h(jj).Color = [0 0 1 0.3]; % fs
                elseif jj>6
                    h(jj).Color = [0 0 0 0.3]; % fns
                end
            end
        else
            ax.XTick=[1:14];
            set(gca,'fontsize',14,'xticklabel',[presynPops(3) presynPops(1) presynPops(2) presynPops(4) " "])
            h = findobj(gca,'Tag','Box');
            for jj=1:length(h)
                if jj<5
                    h(jj).Color = [0.008, 0.424, 0.271 0.3]; % ht
                elseif jj>4 && jj<8
                    h(jj).Color = [0 0 1 0.3]; % fs
                elseif jj>7
                    h(jj).Color = [0 0 0 0.3]; % fns
                end
            end
        end
        title({'Total synaptic input (Full excl.'; ['settling time) per ' char(cellToInspect) ' cell']})
    end
    
    cb = table({sim_tag;sim_tag;sim_tag},{'FNS';'FS';'HT'},...
        {ec_fns;ec_fs;ec_ht},...
        {fns_num;fs_num;ht_num},...
        {faulty_not_spiking;faulty_spiking;healthy},...
        {presyn_ExcPrevLayer.faulty_not_spiking;presyn_ExcPrevLayer.faulty_spiking;presyn_ExcPrevLayer.healthy},...
        {presyn_ExcSameLayer.faulty_not_spiking;presyn_ExcSameLayer.faulty_spiking;presyn_ExcSameLayer.healthy},...
        {presyn_InhSameLayer.faulty_not_spiking;presyn_InhSameLayer.faulty_spiking;presyn_InhSameLayer.healthy},...
        {presyn_InhPrevLayer.faulty_not_spiking;presyn_InhPrevLayer.faulty_spiking;presyn_InhPrevLayer.healthy},...
        {faulty_ns_0presynInh1;faulty_s_0presynInh1;healthy_0presynInh1},...
        {faulty_ns_0presynInh2;faulty_s_0presynInh2;healthy_0presynInh2},...
        {presyn_ExcPrevLayer.faulty_not_spiking(presyn_InhSameLayer.faulty_not_spiking==0);presyn_ExcPrevLayer.faulty_spiking(presyn_InhSameLayer.faulty_spiking==0);presyn_ExcPrevLayer.healthy(presyn_InhSameLayer.healthy==0)},...
        {presyn_ExcSameLayer.faulty_not_spiking(presyn_InhSameLayer.faulty_not_spiking==0);presyn_ExcSameLayer.faulty_spiking(presyn_InhSameLayer.faulty_spiking==0);presyn_ExcSameLayer.healthy(presyn_InhSameLayer.healthy==0)},...
        {presyn_ExcPrevLayer.faulty_not_spiking(presyn_InhSameLayer.faulty_not_spiking~=0);presyn_ExcPrevLayer.faulty_spiking(presyn_InhSameLayer.faulty_spiking~=0);presyn_ExcPrevLayer.healthy(presyn_InhSameLayer.healthy~=0)},...
        {presyn_ExcSameLayer.faulty_not_spiking(presyn_InhSameLayer.faulty_not_spiking~=0);presyn_ExcSameLayer.faulty_spiking(presyn_InhSameLayer.faulty_spiking~=0);presyn_ExcSameLayer.healthy(presyn_InhSameLayer.healthy~=0)},...
        {presyn_ExcPrevLayer.faulty_not_spiking(presyn_InhPrevLayer.faulty_not_spiking==0);presyn_ExcPrevLayer.faulty_spiking(presyn_InhPrevLayer.faulty_spiking==0);presyn_ExcPrevLayer.healthy(presyn_InhPrevLayer.healthy==0)},...
        {presyn_ExcSameLayer.faulty_not_spiking(presyn_InhPrevLayer.faulty_not_spiking==0);presyn_ExcSameLayer.faulty_spiking(presyn_InhPrevLayer.faulty_spiking==0);presyn_ExcSameLayer.healthy(presyn_InhPrevLayer.healthy==0)},...
        {presyn_ExcPrevLayer.faulty_not_spiking(presyn_InhPrevLayer.faulty_not_spiking~=0);presyn_ExcPrevLayer.faulty_spiking(presyn_InhPrevLayer.faulty_spiking~=0);presyn_ExcPrevLayer.healthy(presyn_InhPrevLayer.healthy~=0)},...
        {presyn_ExcSameLayer.faulty_not_spiking(presyn_InhPrevLayer.faulty_not_spiking~=0);presyn_ExcSameLayer.faulty_spiking(presyn_InhPrevLayer.faulty_spiking~=0);presyn_ExcSameLayer.healthy(presyn_InhPrevLayer.healthy~=0)},...
        {tsyn_fns_events_ExcPrevLayer;tsyn_fs_events_ExcPrevLayer;tsyn_ht_events_ExcPrevLayer},...
        {tsyn_fns_events_ExcSameLayer;tsyn_fs_events_ExcSameLayer;tsyn_ht_events_ExcSameLayer},...
        {tsyn_fns_events_InhSameLayer;tsyn_fs_events_InhSameLayer;tsyn_ht_events_InhSameLayer},...
        {tsyn_fns_events_InhPrevLayer;tsyn_fs_events_InhPrevLayer;tsyn_ht_events_InhPrevLayer},...
        {tsyn_fns_full_ExcPrevLayer;tsyn_fs_full_ExcPrevLayer;tsyn_ht_full_ExcPrevLayer},...
        {tsyn_fns_full_ExcSameLayer;tsyn_fs_full_ExcSameLayer;tsyn_ht_full_ExcSameLayer},...
        {tsyn_fns_full_InhSameLayer;tsyn_fs_full_InhSameLayer;tsyn_ht_full_InhSameLayer},...
        {tsyn_fns_full_InhPrevLayer;tsyn_fs_full_InhPrevLayer;tsyn_ht_full_InhPrevLayer},...
        {mem_pot.faulty_not_spiking;mem_pot.faulty_spiking;mem_pot.healthy},...
        {first_spikes.faulty_not_spiking;first_spikes.faulty_spiking;first_spikes.healthy},...
        {STA.faulty_not_spiking;STA.faulty_spiking;STA.healthy},...
        'variablenames',{'sim_tag','behaviour_type','example_cell','cell_count','cell_idx','ExcPrevLayer_conns','ExcSameLayer_conns','InhSameLayer_conns','InhPrevLayer_conns',...
        'noInhSame','noInhPrev','ExcPrevLayer_conns_in_noInhSame','ExcSameLayer_conns_in_noInhSame', 'ExcPrevLayer_conns_in_withInhSame','ExcSameLayer_conns_in_withInhSame',...
        'ExcPrevLayer_conns_in_noInhPrev','ExcSameLayer_conns_in_noInhPrev', 'ExcPrevLayer_conns_in_withInhPrev','ExcSameLayer_conns_in_withInhPrev',...
        'totalSyn_ExcPrevLayer_events','totalSyn_ExcSameLayer_events','totalSyn_InhSameLayer_events','totalSyn_InhPrevLayer_events',...
        'totalSyn_ExcPrevLayer_fullExclSetlling','totalSyn_ExcSameLayer_fullExclSetlling','totalSyn_InhSameLayer_fullExclSetlling','totalSyn_InhPrevLayer_fullExclSetlling',...
        'Vm', 'first_spike_latency', 'STA_vm'});
end

function r = IntegrateEventsSeparately_PerNeuron(neuronData, presynName)
    r = [];
    for evIdx = 2:length(neuronData.(presynName).timeWindows)
        r(evIdx-1) = trapz(neuronData.(presynName).timeWindows{1, evIdx});
    end
    r = sum(r);
end