function results = getSTA_ETA(type,presyn,postsyn,con_mat,bin_size,event_info,win_size,tstart, tstop, isLIF, preParams, postParams)

    if nargin<10
        isLIF = 0;
    end
    
    switch type
        case 'ETA'
            if ~iscell(event_info) && size(event_info,2)~=2 && any(cell2mat(cellfun(@(x) isempty(x), event_info, 'uniformoutput',false)))
                disp('Specify the VPM event detection strategy (threshold/timestamp) and the event threshold/timestamp.')
            else
                strategy = cell2mat(event_info(1));
                threshold_tag = cell2mat(event_info(2));
            end
        case 'STA'
    end
    % get spike times
    tstep = 0.01;
    switch type
        case 'STA'
        % STA (spike-triggered average)
        if ~isLIF
            preSynActivity = getSpikeTimes(presyn, [tstart:tstep:tstop]);
            postSynActivity = getSpikeTimes(postsyn, [tstart:tstep:tstop]);
        else
            preSynActivity = arrayfun(@(x) getLIFspikes(presyn(:,x),[0:tstep:tstop],tstart,preParams(1),preParams(2),0), 1:size(presyn,2),'uniformoutput',false);
            postSynActivity = arrayfun(@(x) getLIFspikes(postsyn(:,x),[0:tstep:tstop],tstart,postParams(1),postParams(2),0), 1:size(postsyn,2),'uniformoutput',false);
        end
        
        % discard spkes before x ms
        x_ms = 100;
        preSynActivity = cellfun(@(x) x(x>x_ms), preSynActivity, 'uniformoutput', false);
        postSynActivity = cellfun(@(x) x(x>x_ms), postSynActivity, 'uniformoutput', false);
        
        % go through all presynaptic neurons, get each spike and look if it triggered a spike in its postsynaptic partner
        results = [];
        for presynNeuron = 1:length(preSynActivity)
            presynSpikeTimes = preSynActivity{1, presynNeuron};  
            if ~isempty(presynSpikeTimes)
                synPartners = find(con_mat(presynNeuron,:)>0);
                synPartners_spikeTimes = postSynActivity(synPartners);
                postsynSpikeTimeLatency = 3; % ms
                % loop through presynaptic spikes
                for preSpike = 1:length(presynSpikeTimes)
                    spiketime = presynSpikeTimes(preSpike);
                    spiketime_latency = spiketime + postsynSpikeTimeLatency;
                    spiketime_preActivation = spiketime - postsynSpikeTimeLatency;
                    postsynResponses_postAP = cellfun(@(x) numel(x(x>spiketime & x<=spiketime_latency)),synPartners_spikeTimes, 'uniformoutput', false);
                    postsynResponses_preAP = cellfun(@(x) numel(x(x>=spiketime_preActivation & x<spiketime)),synPartners_spikeTimes, 'uniformoutput', false);
                    
                    results.presynSpikes(presynNeuron,preSpike) = spiketime;
                    results.postsynSpikes_preAP{presynNeuron,preSpike} = postsynResponses_preAP;
                    results.postsynSpikes_postAP{presynNeuron,preSpike} = postsynResponses_postAP;
                end
            end
        end
        
%%%%%%%%%%%%%%%%%%
        case 'ETA'
        % ETA (event-triggered average)
        VPMactivity = getSpikeTimes(presyn, [tstart:tstep:tstop]);
        if ~isLIF
            CTXactivity = getSpikeTimes(postsyn, [tstart:tstep:tstop]);
        else
            CTXactivity = arrayfun(@(x) getLIFspikes(postsyn(:,x),[0:tstep:tstop],tstart,postParams(1),postParams(2),0), 1:size(postsyn,2),'uniformoutput',false);
        end

        % discard spkes before x ms
        x_ms = 100;
        VPMactivity = cellfun(@(x) x(x>x_ms), VPMactivity, 'uniformoutput', false);
        CTXactivity = cellfun(@(x) x(x>x_ms), CTXactivity, 'uniformoutput', false);
        
        % find VPM events:
        % #1 get vpm psth
        [vpm_psth, vpm_edges] = plotPsthAndRaster('psth',VPMactivity,bin_size,tstop,0); % edges are given in time (ms) with timestep = binSize
        % #2 get peaks in psth
        switch strategy
            case 'threshold'
                [~, vpm_events_location] = findpeaks(vpm_psth,'minpeakheight',threshold_tag,'minpeakdistance',round(win_size/bin_size),'sortstr','descend');
                % adjust events (shift peak if any of the 3 preceeding bins have >spikes_thresh spikes
                %! BEWARE: IF YOU DO THIS, CORTICAL ACTIVITY WILL BE WIDELY
                % DISPERSED BECAUSE THALAMIC EVENTS ARE NOT ALIGNED TO
                % PEAKS BUT TO THE FIRST THRESHOLD CROSSING
                for ev = 1:length(vpm_events_location)
                    if any(vpm_psth([vpm_events_location(ev)-3:vpm_events_location(ev)-1])>threshold_tag)
                        bin_shift_idx = find(vpm_psth([vpm_events_location(ev)-3:vpm_events_location(ev)])>threshold_tag);
                        bin_shift = abs(diff([bin_shift_idx(1) bin_shift_idx(end)]));
                        vpm_events_location(ev) = vpm_events_location(ev)-bin_shift;
                    end
                end; clear ev
                % #3 get time of vpm peaks/events
                time_vec_vpm = vpm_edges(vpm_events_location); 
            case 'peak'
                [~, vpm_events_location] = findpeaks(vpm_psth,'minpeakheight',threshold_tag,'minpeakdistance',round(win_size/bin_size),'sortstr','descend');
                % #3 get time of vpm peaks/events
                time_vec_vpm = vpm_edges(vpm_events_location);                 
            case 'timestamp'
                time_vec_vpm = threshold_tag;
        end
        
        % #4 get ctx psth
        ctx_psth = plotPsthAndRaster('psth',CTXactivity,bin_size,tstop,0);
        % #5 get pre- and post-peak boundaries
        time_pre = 50; % ms pre peak
        time_post = 100; % ms post peak
        pre = round(time_pre/bin_size);
        post = round(time_post/bin_size);
        % #6 loop through vpm events and get ctx psth activity
        results = [];
        results.BinSize = bin_size;
        for ev = 1:length(time_vec_vpm)
            if time_vec_vpm(ev)+time_post > vpm_edges(end)
                last_idx = length(vpm_edges);
                pad_back = zeros(1,abs(diff([time_vec_vpm(ev)+time_post  vpm_edges(end)]))/bin_size);
            else
                last_idx = nearest(vpm_edges, time_vec_vpm(ev))+post;
                pad_back = [];
            end
            if time_vec_vpm(ev)-time_pre < vpm_edges(1)
                first_idx = 1;
                pad_front = zeros(1,abs(diff([time_vec_vpm(ev)-time_post  vpm_edges(1)]))/bin_size);
            else
                first_idx = nearest(vpm_edges, time_vec_vpm(ev)-pre);
                pad_front = [];
            end
            
            results.vpm_events{ev} = [pad_front vpm_psth(first_idx:last_idx) pad_back];
            results.ctx_events{ev} = [pad_front ctx_psth(first_idx:last_idx) pad_back];

            time = [time_vec_vpm(ev)-time_pre time_vec_vpm(ev) time_vec_vpm(ev)+time_post]; % [pre_bound_time peak_time post_bound_time]
            results.event_times{ev} = time;
            
            % get how many cells participate in this event
            EventStartTime = vpm_edges(first_idx);
            EventEndTime = vpm_edges(last_idx);

            CTXactivity_ThisEvent = cellfun(@(x) x(x>EventStartTime & x<EventEndTime), CTXactivity, 'uniformoutput', false);
            VPMactivity_ThisEvent = cellfun(@(x) x(x>EventStartTime & x<EventEndTime), VPMactivity, 'uniformoutput', false);

            results.responsiveCellsNum_vpm{ev} = length(find(cell2mat(cellfun(@(x) ~isempty(x), VPMactivity_ThisEvent, 'uniformoutput', false))>0));
            results.responsiveCellsNum_ctx{ev} = length(find(cell2mat(cellfun(@(x) ~isempty(x), CTXactivity_ThisEvent, 'uniformoutput', false))>0));
            
%             % spike rate
%                 % FR pre event
%             results.spike_rate_pre_vpm{ev} = cell2mat(cellfun(@(x) numel(x(x>=time(1) & x<time(2)))/(diff(time([1 2])/1000)), VPMactivity, 'uniformoutput', false));
%             results.spike_rate_pre_ctx{ev} = cell2mat(cellfun(@(x) numel(x(x>=time(1) & x<time(2)))/(diff(time([1 2])/1000)), CTXactivity, 'uniformoutput', false));
%                 % FR post event
%             results.spike_rate_post_vpm{ev} = cell2mat(cellfun(@(x) numel(x(x>=time(2) & x<time(3)))/(diff(time([2 3])/1000)), VPMactivity, 'uniformoutput', false));
%             results.spike_rate_post_ctx{ev} = cell2mat(cellfun(@(x) numel(x(x>=time(2) & x<time(3)))/(diff(time([2 3])/1000)), CTXactivity, 'uniformoutput', false));  
%             
%             % inverse isi
%                 % FR pre event
%             results.inst_fr_pre_vpm{ev} = cell2mat(cellfun(@(x) 1./(diff(x(x>=time(1) & x<time(2)))/1000), VPMactivity, 'uniformoutput', false));
%             results.inst_fr_pre_ctx{ev} = cell2mat(cellfun(@(x) 1./(diff(x(x>=time(1) & x<time(2)))/1000), CTXactivity, 'uniformoutput', false));
%                 % FR post event
%             results.inst_fr_post_vpm{ev} = cell2mat(cellfun(@(x) 1./(diff(x(x>=time(2) & x<time(3)))/1000), VPMactivity, 'uniformoutput', false));
%             results.inst_fr_post_ctx{ev} = cell2mat(cellfun(@(x) 1./(diff(x(x>=time(2) & x<time(3)))/1000), CTXactivity, 'uniformoutput', false));               
        

            % latency
            ltc_first = []; ltc_second = [];
            for nrn = 1:size(CTXactivity,2)
                first_spike_time = CTXactivity{1, nrn}(find(CTXactivity{1, nrn}>time_vec_vpm(ev),1,'first'));
                second_spike_time = CTXactivity{1, nrn}(find(CTXactivity{1, nrn}>time_vec_vpm(ev),2,'first'));
                
                if isempty(first_spike_time)
                    first_spike_time = NaN;
                    second_spike_time = NaN;
                end
                if length(second_spike_time)<=1
                    second_spike_time = NaN;
                else
                    second_spike_time = second_spike_time(2);
                end
                                
                ltc_first(nrn,1) = first_spike_time - time_vec_vpm(ev);
                ltc_second(nrn,1) = second_spike_time - time_vec_vpm(ev);
            end
            results.first_spike_latency{ev} = ltc_first;
            results.second_spike_latency{ev} = ltc_second;
        end
    end
end