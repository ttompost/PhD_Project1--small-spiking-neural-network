function spike_times = getSpikeTimes(s, time)
%% this is a method from dsCalcFR.m .. not really because it did not work so I modified it
%     time = [100:0.01:1100]; % 1 s simulation with 100ms discarded in the beginning 
    threshold = -18;
    ncells = size(s,2);
    spike_times=cell(1,ncells);
    for ii=1:ncells
        % get spikes in this cell
        spike_inds_left=1+find((s(2:end,ii)>=threshold & s(1:end-1,ii)<threshold));
        spike_inds_right = spike_inds_left+(0.7/0.01);
        spike_inds = [];
        for jj = 1:length(spike_inds_left)
            try
                mx = max(s(spike_inds_left(jj):spike_inds_right(jj),ii));
                spike_inds(jj) = find(s(spike_inds_left(jj):spike_inds_right(jj),ii) == mx);
            catch ME
                spike_inds(jj) = find(s(spike_inds_left(jj),ii));
            end
            
        end
        spike_inds = spike_inds+spike_inds_left'-1;
        spike_times{ii}=time(spike_inds);
    end
end