function results = plotAPs(AP_struct, stimAmps, stimOn, stimDur, simOff, col, whichToPlot)

    stimTime = stimDur/1000; % in sec
    plot_tag = 0;
    
    if isempty(col) && ~isempty(whichToPlot)
        colors = flipud(colormap(winter(length(whichToPlot)+3)));
    elseif isempty(col) && isempty(whichToPlot)
        colors = flipud(colormap(hot(length(stimAmps)+3)));
    end
    if isempty(whichToPlot)
        whichToPlot = 1:size(AP_struct,2);
    end
    lgd = [];
    results = table([],[],[],[],[],[],[],[],[],'variablenames', {'spike_times_full';'spike_times_range';'ISI_mean_full'; 'ISI_std_full'; 'ISI_mean_range'; 'ISI_std_range'; 'FR_full'; 'FR_range'; 'Iapp'});
    for ii = 1:size(AP_struct,2)
        if ismember(ii, whichToPlot)
            if exist('colors','var')
                plot(0:0.01:simOff,AP_struct(ii).(AP_struct(ii).labels{1, 1})(:,1), 'color', colors(ii,:), 'linewidth', 1)   
            else
                plot(0:0.01:simOff,AP_struct(ii).(AP_struct(ii).labels{1, 1})(:,1), 'color', col, 'linewidth', 1) 
            end
            lgd(end+1) = stimAmps(ii);
            plot_tag = 1;
        end
        
        spikes = cell2mat(AP_struct(ii).([AP_struct(ii).labels{1, 1} '_spike_times']));
        spikes_range = spikes(spikes>stimOn & spikes<stimOn+stimDur);
        
        isi_mean_full = mean(diff(spikes));
        isi_std_full = std(diff(spikes));
        
        isi_mean_range = mean(diff(spikes_range));
        isi_std_range = std(diff(spikes_range));
        
        fr_calc_range = numel(spikes_range)/stimTime; % only spiked in stimulation time
        fr_calc_full = numel(spikes)/simOff; % all spikes in full time
        
        results = [table({spikes},{spikes_range},{isi_mean_full},{isi_std_full},{isi_mean_range},{isi_std_range},{fr_calc_full},{fr_calc_range},{stimAmps(ii)}, 'variablenames', ...
            {'spike_times_full';'spike_times_range';'ISI_mean_full'; 'ISI_std_full'; 'ISI_mean_range'; 'ISI_std_range'; 'FR_full'; 'FR_range'; 'Iapp'}); results];
    end
   
    if plot_tag
        l = legend(string(lgd), 'location','eastoutside','fontsize', 10);
        title(l, {'Current'; 'amplitude (pA)'})
        ylim([-100 50])
        xlim([0 simOff])
        set(gca,'fontsize', 14)
        xlabel('Time (ms)')
        ylabel('Vm (mV)')
    end
end