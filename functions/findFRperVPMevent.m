function events = findFRperVPMevent(VPMactivity,CTXactivity, Bin, Spikes, Win, Lay)
    % VPMactivity = cell array of spike timings per VPM neuron
    % CTXactivity = cell array of spike timings per cortical neuron
    % Bin = bin size for PSTH and event detection
    % Spikes = number of spikes per bin to define an event
    % Win = window size for FR calculation
    % Lay = layer number

    % find psth and peaks
    [vpm_psth, vpm_edges] = plotPsthAndRaster('psth',VPMactivity,Bin,0,[]);
    [~, vpm_events_location] = findpeaks(vpm_psth,'minpeakheight',Spikes,'minpeakdistance',round(Win/Bin),'sortstr','descend');
    time_vec_vpm = vpm_edges(vpm_events_location-1); % shift some bins before
    if Lay == 4
        time_vec_ctx = time_vec_vpm + 1; % l4 response is delayed for 1 bin (of size 2ms); however, that is true for l4i. l4e are delayed for 1ms, therefore, their bins are aligned with vpm...
    elseif Lay == 3
        time_vec_ctx = time_vec_vpm + 2; % 2 bins * 2ms delay
    elseif Lay == 2
        time_vec_ctx = time_vec_vpm + 3; % 3 bins * 2ms delay
    else
        disp('Specify layer. Terminating funciton..')
        return
    end

    % find APs per event
    for event = 1:length(time_vec_ctx)
        aps = [];
        for ctx_neuron = 1:size(CTXactivity,2)
            aps(ctx_neuron) = numel(find(CTXactivity{1, ctx_neuron} > time_vec_ctx(event) & CTXactivity{1, ctx_neuron} < time_vec_ctx(event) + Win));
        end
        events{event,1} = [time_vec_ctx(event) time_vec_ctx(event) + Win];
        events{event,2} = aps;
        
        aps = [];
        for vpm_neuron = 1:size(VPMactivity,2)
            aps(vpm_neuron) = numel(find(VPMactivity{1, vpm_neuron} > time_vec_vpm(event) & VPMactivity{1, vpm_neuron} < time_vec_vpm(event) + Win));
        end
        events{event,3} = [time_vec_vpm(event) time_vec_vpm(event) + Win];
        events{event,4} = aps;
    end
end