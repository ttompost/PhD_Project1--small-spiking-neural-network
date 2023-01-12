function p = getPeriodicityFromAutox(type, sig, minheight, mindist, T, smth)
    % sig has to be in raw datapoints
    % mindist has to be in milliseconds
    % T is sampling period in ms
    
    % mean peak height will be used as min peak prominence
    warning('off')
    half_sig = sig(round(length(sig)/2) : end); % get positive half of autox signal
    switch type
        case 'vm'
            lgth = floor(length(sig)*T)/2;
            tvec=-lgth:T:lgth;
            if nargin == 6
                half_sig = movmean(half_sig,smth/T);
            end
        case {'psth';'psth_single'}
            tvec=-(length(sig)-1):T:(length(sig)-1);
    end
    half_tvec = tvec(find(tvec==0):end);

% % PREVIOUS APPROACH
%     [pks, locs] = findpeaks(half_sig,'minpeakprominence',minheight,'minpeakdist',mindist/T,'npeaks', 3, 'minpeakheight', minheight, 'sortstr', 'none'); % find first 3 peaks
%     if ~isempty(locs)
%         locs_ms = [0 locs]*T; % add first peak at 0 and turn location indices into miliseconds
%         if std(diff(locs_ms)) > 10 % ms
%             if std(diff(locs_ms(1:end-1))) > 5 % remove last peak to see if it helps std
%                 [pks, locs] = findpeaks(half_sig,'minpeakprominence',minheight,'minpeakdist',(mindist/T)*4,'npeaks', 3, 'minpeakheight', minheight, 'sortstr', 'desc');
%                 locs_ms = [0 sort(locs)]*T; % add first peak at 0 and turn location indices into miliseconds
%                 mean_difflocs = mean(diff(locs_ms)); % get difference in locations, i.e. how often the peak occurs, and then get its mean
%                 p = 1/(mean_difflocs/1000);
%                  if std(diff(locs_ms)) > 10
%                     p=Inf; % if detected peaks have big std, then there is something wrong with detection (or there is no periodicity)
%                  end
%             else
%                 locs_ms = [0 locs(1:end-1)]*T; % add first peak at 0 and turn location indices into miliseconds
%                 mean_difflocs = mean(diff(locs_ms)); % get difference in locations, i.e. how often the peak occurs, and then get its mean
%                 p = 1/(mean_difflocs/1000);
%             end
%         else
%             mean_difflocs = mean(diff(locs_ms)); % get difference in locations, i.e. how often the peak occurs, and then get its mean
%             p = 1/(mean_difflocs/1000);
%         end
%     else
%         p=Inf; % signal does not oscillate, infinite periodicity
%     end

% % NEW APPROACH
    [~, locs1, wdth1] = findpeaks(sig,'minpeakprominence',minheight,'npeaks', 1, 'minpeakheight', minheight, 'sortstr', 'desc', 'WidthReference','halfprom');
    half_WidthZero = round(ceil(wdth1)/2);
    PeakAtZero = sig(locs1-half_WidthZero:locs1+half_WidthZero);
    time1 = tvec(locs1); % in ms

    [~, locs2, wdth2] = findpeaks(half_sig,'minpeakprominence',minheight,'npeaks', 1, 'sortstr', 'none', 'WidthReference','halfprom');
    if ~isempty(wdth2) && locs2*T < 200
        half_WidthFirst = round(ceil(wdth2)/2);
        if (locs2-half_WidthFirst) < 0 || (locs2+half_WidthFirst) > length(half_sig)
            [~, locs2, wdth2] = findpeaks(half_sig,'minpeakprominence',minheight,'npeaks', 1, 'sortstr', 'none', 'WidthReference','halfheight');
            half_WidthFirst = round(ceil(wdth2)/2);
            if half_WidthFirst*T > 35
                FirstPeak = half_sig(locs2-half_WidthFirst:end);
            else
                FirstPeak = half_sig(locs2-half_WidthFirst:locs2+half_WidthFirst);
            end
        else
            FirstPeak = half_sig(locs2-half_WidthFirst:locs2+half_WidthFirst);
        end
    else
        time2 = 500;
    end
    
    try
        flank_left = half_sig(locs2-half_WidthFirst-length(FirstPeak):locs2-half_WidthFirst);
        if half_WidthFirst*T < 35
            flank_right = half_sig(locs2+half_WidthFirst:locs2+half_WidthFirst+length(FirstPeak));
        end
        time2 = half_tvec(locs2); % in ms
    catch ME
        time2 = 500; % or whatever is bigger than the number in following if-cond
    end
    
    if time2 < 200 % not interested in peaks after 200ms
        power_PAZ = bandpower(PeakAtZero); % not sure what to do with this yet
        power_FP = bandpower(FirstPeak);
        power_FL = bandpower(flank_left);
        if half_WidthFirst*T < 35
            power_FR = bandpower(flank_right);
        else
            power_FR = NaN;
        end
        std_tail = std(half_sig(round(length(half_sig)/2):end));
        
        snr_left = power_FP/power_FL;
        snr_right = power_FP/power_FR;
        if all([snr_left snr_right] > std_tail/2)
            p = 1/(diff([time1 time2])/1000);
        elseif any([snr_left snr_right] > std_tail)
            p = 1/(diff([time1 time2])/1000);
        else
            p = Inf;
        end
    else 
        p = Inf;
    end
    
    warning('on')
end