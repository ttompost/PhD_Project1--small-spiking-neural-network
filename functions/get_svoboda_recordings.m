function [whiskermat,binsize_whisker, whisker_metadata] = get_svoboda_recordings(loadfolder, animal, sessionvec, dataname, window)
% this function gets whisker angle in the 1st dimension, and whisker curvature
% in the 2nd dimention

% window has to have:
%       window.start = ''
%       window.window = [start, stop]

if nargin < 5
    window.window = [];
    window.start = [];
elseif nargin == 5
    if isnumeric(window)
        window.window = [];
        window.start = [];
    end
end

whisker_metadata = table([],[],[],[],[],[],[],[],[],[],'variablenames',{'AnimalName','SessionDate','TrialNum','PoleLocation_mm','BehavioralResponse','WhiskerRecording','TouchDetected','WhiskerAngle_deg','WhiskerCurve_mm-1','time'});

Nsession = length(sessionvec);
for ns = 1:Nsession
    filename = [loadfolder char(animal) '_' sessionvec{ns} '_' dataname '.mat'];
    disp(['Loading file ' filename])

    load(filename) % any file from https://crcns.org/data-sets/ssc/ssc-2

    % check if multiple whiskers were recorded (how can multiple whiskers be recorded in a single-whisker experiment?)
    [nwhiskertrace, ~] = size(s.timeSeriesArrayHash.value{1}.valueMatrix);
    if nwhiskertrace>2
        disp('More than 2 whisker sessions detected:')
        for nw = 1:nwhiskertrace
            disp(s.timeSeriesArrayHash.value{1}.idStrs{nw})
        end
%         whiskertracen = input('Which traces should be used? (give 2d array)');

        % if i understand correctly, when there are multiple whisker
        % recordings, then they are sorted as: all angle traces first, and
        % then all curvature traces.. so i will take the first recording of
        % each
        ang_idx = find(cell2mat(cellfun(@(x) contains(char(x), 'Angle'),s.timeSeriesArrayHash.value{1, 1}.idStrs, 'uniformoutput',false)),1,'first')+1;
        curv_idx = find(cell2mat(cellfun(@(x) contains(char(x), 'Curvature'),s.timeSeriesArrayHash.value{1, 1}.idStrs, 'uniformoutput',false)),1,'first')+1;
        whiskertracen = [ang_idx,curv_idx];
        
        disp(['I am using the second session. Angle, ' num2str(ang_idx) '; and curvature, ' num2str(curv_idx) '.'])
    else
        whiskertracen = [1,2];
    end

    % number of trials for this volume
    volume = 2; % volume matters if you want to consider calcium traces, but for whisker traces it does not matter. i think.
    trialvec = unique(s.timeSeriesArrayHash.value{volume}.trial);
    Ntrial_temp = length(trialvec);
    binsize_whisker = s.timeSeriesArrayHash.value{1}.time(2)-s.timeSeriesArrayHash.value{1}.time(1);

    if strcmp(window.start, 'pole in reach')
        Ntime_w = 1000;
    elseif (~isempty(window.window) && ~strcmp(window.start, 'pole in reach'))
        whiskertime = window.window(1):binsize_whisker:window.window(end);
        Ntime_w = length(whiskertime);
    else
        error('Please give an appropriate time window')
    end

    whiskermat_temp = NaN*ones(2,Ntime_w,Ntrial_temp);

    true_trial_idx = 0;
    skip_trial = 0;
    disp(['There are ' num2str(Ntrial_temp) ' trials in this session.'])
    for nt = 1:Ntrial_temp
        whiskertrialvec = find(s.timeSeriesArrayHash.value{1}.trial == trialvec(nt));
        luminestrialvec = find(s.timeSeriesArrayHash.value{volume}.trial == trialvec(nt));
        
        poleLocation_Zaber = s.trialPropertiesHash.value{1, 3}(trialvec(nt));
        poleLocation_mm = poleLocation_Zaber/10000; % 1 mm = 10000 zaber motor steps
        
        % is the trial a hit
        behavResponse_idx = find(s.trialTypeMat(:,trialvec(nt))); 
        behavResponse = cell2mat(s.trialTypeStr(behavResponse_idx));
        
        if contains(behavResponse,'Hit')
            is_a_hit = 1;
        else
            is_a_hit = 0;
        end
        

        % extract whiskking data
        if isempty(whiskertrialvec) || any(isnan(whiskertrialvec))
            disp(['No whisker data, skipping trial with id ' num2str(trialvec(nt))])
            
            whisker_rec = 0;
            touch_rec = 0;
%         elseif ~is_a_hit % I will first ignore this condition, it turns out that not a lot of them are a hit
%             disp(['Trial was not a hit, skipping trial with id ' num2str(trialvec(nt))])

        else
            whisker_rec = 1;
            
            whiskertime_thistrial = s.timeSeriesArrayHash.value{1}.time(whiskertrialvec);
            luminestime_thistrial = s.timeSeriesArrayHash.value{volume}.time(luminestrialvec);

            whiskertrace = s.timeSeriesArrayHash.value{1}.valueMatrix(whiskertracen,whiskertrialvec);

            tstart = max(whiskertime_thistrial(1), luminestime_thistrial(1));
            tend = min(whiskertime_thistrial(end), luminestime_thistrial(end));

            if isempty(window.start)
               % nothing needed, use calculated tstart and tend;
               touch_rec = 1;
            elseif strcmp(window.start, 'pole in reach')
                if isfield(window, 'window')
                    disp('Using times pole in reach; ignoring given window')
                end
                tpole  = s.eventSeriesArrayHash.value{1}.eventTimes(s.eventSeriesArrayHash.value{1}.eventTrials==trialvec(nt));
                tstart = tpole(1);
                tend   = tpole(2);
                
                touch_rec = 1;
            elseif strcmp(window.start, 'first touch')
                ttouchpro = s.eventSeriesArrayHash.value{2}.eventTimes{1}(s.eventSeriesArrayHash.value{2}.eventTrials{1}==trialvec(nt));
                ttouchre  = s.eventSeriesArrayHash.value{2}.eventTimes{2}(s.eventSeriesArrayHash.value{2}.eventTrials{2}==trialvec(nt));
                if ~isempty(ttouchpro) && ~isempty(ttouchre)
                    ttouchpro = ttouchpro(1);
                    ttouchre  = ttouchre(1);
                    ttouch = min(ttouchpro, ttouchre);
                    
                    touch_rec = 1;
                elseif ~isempty(ttouchpro) && isempty(ttouchre)
                    ttouch = ttouchpro(1);
                    
                    touch_rec = 1;
                elseif isempty(ttouchpro) && ~isempty(ttouchre)
                    ttouch = ttouchre(1);
                    
                    touch_rec = 1;
                else
%                     disp(['No touch data, skipping trial with id ' num2str(trialvec(nt))])
                    ttouch = tstart;
                    skip_trial = skip_trial + 1;
                    
                    touch_rec = 0;
                end                    
                tstart = ttouch+window.window(1);
                tend   = ttouch+window.window(2);
            elseif strcmp(window.start, 'first')
                tstart = tstart+window.window(1);
                tend   = tstart+window.window(2);
                
                touch_rec = 1;
            end

            if strcmp(window.start, 'pole in reach')
                % variable length
                [~, nstart] = min(abs(whiskertime_thistrial-tstart));
                [~, nend] = min(abs(whiskertime_thistrial-tend));
                whiskertrace = whiskertrace(:,nstart:nend);
                
                % register a true trial now
                true_trial_idx = true_trial_idx+1;
            else
                % fixed length
                if ((tstart>whiskertime_thistrial(1)) && (tend<whiskertime_thistrial(end)))
                    % window fits in trial
                    [~, nstart] = min(abs(whiskertime_thistrial-tstart));
                    whiskertrace = whiskertrace(:,nstart:nstart+Ntime_w-1);
                    
                    whiskertrace(1,:) = fillmissing(whiskertrace(1,:), 'linear');
                    whiskertrace(2,:) = fillmissing(whiskertrace(2,:), 'linear');
                    
                    % register a true trial now
                    true_trial_idx = true_trial_idx+1;
                elseif (tstart<whiskertime_thistrial(1) && ~strcmp(window.start, 'pole in reach'))
                    % no recording in beginning of window, add nan 
                    [~, nend] = min(abs(whiskertime_thistrial-tend));
                    whiskertrace_temp = whiskertrace(:,1:nend);
                    lw = length(whiskertrace_temp(1,:));
                    whiskertrace = [nan*ones(2,Ntime_w-lw), whiskertrace_temp];
                    
                    notNaN1 = whiskertrace(1,find(~isnan(whiskertrace(1,:)), 1, 'first'));
                    notNaN2 = whiskertrace(2,find(~isnan(whiskertrace(2,:)), 1, 'first'));
                    
                    whiskertrace(1,1) = notNaN1;
                    whiskertrace(2,1) = notNaN2;
                    
                    whiskertrace(1,:) = fillmissing(whiskertrace(1,:), 'linear');
                    whiskertrace(2,:) = fillmissing(whiskertrace(2,:), 'linear');
                    
                    % register a true trial now
                    true_trial_idx = true_trial_idx+1;
                elseif (tend>whiskertime_thistrial(end) && ~strcmp(window.start, 'pole in reach'))
                    % no recording in end of window, add zeros at the nan
                    [~, nstart] = min(abs(whiskertime_thistrial-tstart));
                    whiskertrace_temp = whiskertrace(:,nstart:end);
                    lw = length(whiskertrace_temp(1,:));
                    whiskertrace = [whiskertrace_temp, nan*ones(2,Ntime_w-lw)];
                    
                    whiskertrace(1,:) = fillmissing(whiskertrace(1,:), 'linear');
                    whiskertrace(2,:) = fillmissing(whiskertrace(2,:), 'linear');
                    
                    % register a true trial now
                    true_trial_idx = true_trial_idx+1;
                else
                    error('Chosen window too large for trials, please chose a smaller window')
                end
            end
            
            % center the signals at zero
            if whiskertrace(1,1) ~= 0 
                whiskertrace(1,:) = whiskertrace(1,:)-whiskertrace(1,1);
            end
            if whiskertrace(2,1) ~= 0
                whiskertrace(2,:) = whiskertrace(2,:)-whiskertrace(2,1);
            end
            
            % standardize signals
            if any(whiskertrace(1,:)>0) || any(whiskertrace(1,:)>0)
                whiskertrace(1,:) = whiskertrace(1,:)/std(whiskertrace(1,:));
            end
            if any(whiskertrace(2,:)>0) || any(whiskertrace(2,:)>0)
                whiskertrace(2,:) = whiskertrace(2,:)/std(whiskertrace(2,:));
            end
            
            if strcmp(window.start, 'pole in reach')
                nwl = length(whiskertrace(1,:));
                if nwl>nwlmax
                    nwlmax = nwl;
                end
                whiskermat_temp(:,1:nwl,nt) = whiskertrace;
            else
                try
                    whiskermat_temp(:,:,nt) = whiskertrace;
                catch
                    keyboard
                end
            end
        end
        whisker_metadata = [whisker_metadata; table(animal,sessionvec(ns),trialvec(nt),poleLocation_mm,{behavResponse},whisker_rec,touch_rec,{whiskertrace(1,:)},{whiskertrace(2,:)},{whiskertime},...
            'variablenames',{'AnimalName','SessionDate','TrialNum','PoleLocation_mm','BehavioralResponse','WhiskerRecording','TouchDetected','WhiskerAngle_deg','WhiskerCurve_mm-1','time'})];
    end
    disp(['I skipped ' num2str([Ntrial_temp-true_trial_idx]) ' trials, and there were ' num2str(skip_trial) ' no-touch trials. They were supposed to be skipped? Not sure.'])
    
    % concatenate with previous sessions
    if exist('whiskermat','var')
        if ~isempty(whiskermat_temp)
            whiskermat = cat(3,whiskermat, whiskermat_temp);
        end
    else
        whiskermat = whiskermat_temp;
    end
    
end

%%%
end
