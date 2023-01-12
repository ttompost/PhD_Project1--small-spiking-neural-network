function ipscAmps = findIpscSize(s)
    
    % maybe these params should be available to the user, but not for now
    dt=0.01; 
    lag = 5; % msec
    peak_influence = 1;
    std_increase = 3; % zscore
    least_amp = 3; % mV
    ipsc_time_window = 0.5; % msec
    same_ipsc = 0.5; % msec
    ipsc_thr = -1.2; % mV 

    ipsc_peaks = threshAlg(s, lag/dt, std_increase, peak_influence);
    ipsc_peaks = find(ipsc_peaks == -1);

    % find consecutive numbers
    k = 1; % to store ipsc amps
    n = 1; % to count the island size
    island = [];% iportion of a signal that crosses zscore boundary, out of which only one datapoint is the ipsc peak
    ipscAmps = [];

    for ii = 1:length(ipsc_peaks)-1
        if s(ipsc_peaks(ii)) < ipsc_thr
            if diff([ipsc_peaks(ii),ipsc_peaks(ii+1)]) < same_ipsc/dt % every peak detected in this window is the same ipsc 
                island([n n+1]) = [ipsc_peaks(ii) ipsc_peaks(ii+1)];
                n = n+1;
            else
                if ~isempty(island) % here i should have an island of consecutive indices to work with
                    i_peak = min(s(island));
                    i_idx = find(s == i_peak);
    
                    if length(i_idx) > 1
                        i_idx = island(ismember(island, i_idx));
                        if length(i_idx) > 1
                            i_idx = i_idx(1);
                        end
                    end
    
                    i_amp = abs(diff([i_peak max(s(i_idx-(ipsc_time_window/dt):i_idx))])); % amplitude of the ipsc
                    if i_amp > least_amp
                        ipscAmps(k,1) = i_amp;
                        ipscAmps(k,2) = i_idx;
                        k = k+1;
                    end
                    island = []; n = 1;
                end
            end
        end
    end
    if isempty(ipscAmps)
        ipscAmps = NaN;
    end
end