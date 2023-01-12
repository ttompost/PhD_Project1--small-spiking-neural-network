function epspAmps = findEpspSize(s, std_increase)
    % provide a membrane potential signal for epsp detection and analysis
    % (for now it only outputs epsp amplitude)
    
    % make sure that the signal is truncated according to the time it takes
    % for the network to settle (e.g., simulation time = 1.1s, truncated
    % part is the initial 0.1s)
    warning('off')
    dt = 0.01;
    smoothing_window = 20; % in ms
    cutAps = s;
    cutAps(cutAps > -55) = -55;
    smooth_s = movmean(cutAps, smoothing_window/dt); % only used to calculate the static event threshold
    peak_distance = 1; % in ms
    event_threshold = -18; % for AP
    peak_prom = 3;
    [pA,pI]=findpeaks(s,'minpeakheight',event_threshold,'MinPeakDistance',peak_distance/dt, 'MinPeakProminence', peak_prom);
    epspAmps = table([],[],[],[],'VariableNames',{'peakIdx','peakAmp','apThreshold','epspAmp'});

    if ~isempty(pA) % if there are action potentials in the signal
        for ap_num = 1:length(pA)
            if pI(ap_num) > 1/dt && pI(ap_num) < 999/dt % ap must not occur in the very first or the very last part because then i am out of array boundaries when taking a window
                bef = ((pI(ap_num) * dt) - (peak_distance/2))/dt; % time before the peak
                aft = ((pI(ap_num) * dt) + (peak_distance/4))/dt; % time after the peak
                ap = s(bef:aft); % action potential
    
                % find dAP/dt
                dAP = diff(ap);
    
                % find d2AP/dt2
                ddAP = diff(dAP);
                Vt = max(ddAP); % spike threshold aka EPSP amplitude
               
                epspAmps = [table(pI(ap_num),pA(ap_num),Vt,NaN(size(Vt)),'VariableNames',{'peakIdx','peakAmp','apThreshold','epspAmp'}); epspAmps];
            end
        end
        % but also check for epsps which failed to cross the threshold
        % std_increase = 2; % best is: 0.5 for I cells; 4 for E cells   % this is  now being given as an input variable
%         peak_distance = 5;
        static_thr = mean(movmedian(smooth_s,40/dt)) + (std_increase * mean(movstd(smooth_s,40/dt))); %static threshold
        [pA_d, pI_d, ~, pProms_d] = findpeaks(s, 'minpeakheight', static_thr, 'minpeakdistance', peak_distance/dt, 'MinPeakProminence', peak_prom);
        epspAmps = [table(pI_d, pA_d, NaN(size(pI_d)), pProms_d,'VariableNames',{'peakIdx','peakAmp','apThreshold','epspAmp'}); epspAmps]; % store additional events

    else % if there are no action potentials in the signal
        std_increase_factor = 4; % for how many STDs should the threshold increase
        peak_distance = 3; % in  ms
        %peak_prom = 0.5;

        static_thr = mean(movmedian(smooth_s,40/dt)) + (std_increase_factor * mean(movstd(smooth_s,40/dt))); %static threshold

        [pA, pI, ~, pProms] = findpeaks(s, 'minpeakheight', static_thr, 'minpeakdistance', peak_distance/dt, 'MinPeakProminence', peak_prom); 
       
        epspAmps = [table(pI, pA, NaN(size(pI)), pProms, 'VariableNames',{'peakIdx','peakAmp','apThreshold','epspAmp'}); epspAmps];
    end


    if size(epspAmps,1) == 0
        disp('no events detected')
    end
    
    % filter for only epsps, not aps
    if any(find(epspAmps.peakAmp > - 55))
        epspAmps(find(epspAmps.peakAmp > - 55),:) = [];
    end

    epspAmps; %  this is just to be able to pause the function 
    warning('on')
end