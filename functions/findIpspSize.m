function ipspAmps = findIpspSize(s)
    
    % maybe these params should be available to the user, but not for now
    dt=0.01; 
    lag = 3.2; % msec
    peak_influence = 1;
    std_increase = 2; % zscore
    least_amp = 4.5; % mV
    ipsp_time_window = 1.2; % msec // although, ipsp rise time is 1 ms, but i increased the window in case ipsp peak is detected later
    same_ipsp = 2; % msec
    ipsp_thr = -63; % mV // if peak is detected above that, it's definitely not an ipsp, but an ap decay in a burst or something

    ipsp_peaks = threshAlg(s, lag/dt, std_increase, peak_influence);
    ipsp_peaks = find(ipsp_peaks == -1);

    % find consecutive numbers
    k = 1; % to store ipsp amps
    n = 1; % to count the island size
    island = [];% portion of a signal that crosses zscore boundary, out of which only one datapoint is the ipsp peak
    ipspAmps = [];

    for ii = 1:length(ipsp_peaks)-1
        if s(ipsp_peaks(ii)) < ipsp_thr
            if diff([ipsp_peaks(ii),ipsp_peaks(ii+1)]) < same_ipsp/dt % every peak detected in this window is the same ipsp 
                island([n n+1]) = [ipsp_peaks(ii) ipsp_peaks(ii+1)];
                n = n+1;
            else
                if ~isempty(island) % here i should have an island of consecutive numbers to work with
                    i_peak = min(s(island));
                    i_idx = find(s == i_peak);
    
                    if length(i_idx) > 1
                        i_idx = island(ismember(island, i_idx));
                        if length(i_idx) > 1
                            i_idx = i_idx(1);
                        end
                    end
    
                    i_amp = abs(diff([i_peak max(s(i_idx-(ipsp_time_window/dt):i_idx))])); % amplitude of the ipsp
                    if i_amp > least_amp
                        ipspAmps(k,1) = i_amp;
                        ipspAmps(k,2) = i_idx;
                        k = k+1;
                    end
                    island = []; n = 1;
                else % if i dont have an island, i take the value of one single point
                    i_peak = s(ipsp_peaks(ii));
                    i_idx = ipsp_peaks(ii);

                    i_amp = abs(diff([i_peak max(s(i_idx-(ipsp_time_window/dt):i_idx))])); % amplitude of the ipsp
                    if i_amp > least_amp
                        ipspAmps(k,1) = i_amp;
                        ipspAmps(k,2) = i_idx;
                        k = k+1;
                    end
                end
            end
        end
    end
    if isempty(ipspAmps)
        ipspAmps = NaN;
    else
        if any(find(ipspAmps(:,1) > 40))
%            p % something not right, pause the function
           ipspAmps(find(ipspAmps(:,1) > 40),:) = [];
        end
    end
end