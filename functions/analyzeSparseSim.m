function AnalysedResults = analyzeSparseSim(Model,lbl)
    %%% 
    % Function description..
    %%%
    
    AnalysedResults = [];
    try
        SeedNum = Model.SeedNum;
        ScalingFactor = [num2str(Model.SparseExpInfo{1, 3}*100) '%'];
        ScaledParameter = Model.SparseExpInfo{1, 2};
        
        GotASeed = 1;
    catch ME
        GotASeed = 0;
    end
    
    if ~GotASeed
        try
            SeedNum = Model.Seeds;
            ScalingFactor = [];
            ScaledParameter = [];
        catch ME
        end
    else
        SeedNum = [];
        ScalingFactor = [];
        ScaledParameter = [];
    end
    
    %% Population names
    PopNamesLong = Model.labels(find(cell2mat(cellfun(@(x) endsWith(x,'_V'), Model.labels, 'uniformoutput', false))));
    PopNamesSplit = cellfun(@(x) split(x, '_'), PopNamesLong, 'uniformoutput', false);
    PopNames = cellfun(@(x) x{1,1}, PopNamesSplit, 'uniformoutput', false);
    
    CTXNamesLong = PopNamesLong(find(cell2mat((cellfun(@(x) startsWith(x, 'L'), PopNamesLong, 'uniformoutput', false)))));
    CTXNames = PopNames(find(cell2mat((cellfun(@(x) startsWith(x, 'L'), PopNames, 'uniformoutput', false)))));
    
    VPMNamesLong = PopNamesLong(find(cell2mat((cellfun(@(x) startsWith(x, 'V'), PopNamesLong, 'uniformoutput', false)))));

    %% Internal Dynasim function to get spike-times
    try
        Model=dsCalcFR(Model);
    catch ME
    end
    CutStartTime = 100; % I do not consider activity that happens before 100ms; the network is settling at t<100ms
    SimulationTime_InMs = Model.time(end)-CutStartTime;
    SimulationTime_InSec = SimulationTime_InMs/1000;
    
    SamplingPeriod = 0.01; % in ms
    CutTimeVector = CutStartTime:SamplingPeriod:Model.time(end);
    
    %% Detect thalamic (VPM) events
    % VPM neuron activity
    CombinedVPM_Spikes = [Model.VPMuncorr_V_spike_times Model.VPMcorr_V_spike_times];
    FilteredVPM_Spikes = cellfun(@(x) x(x>CutStartTime), CombinedVPM_Spikes, 'uniformoutput', false);

    PSTHBinSize = 3; %ms
    [PSTH_VPMPop, Edges_VPMPop] = plotPsthAndRaster('psth',FilteredVPM_Spikes,PSTHBinSize,SimulationTime_InMs,0);

    MinimalEventHeight = mean(PSTH_VPMPop) + (1*std(PSTH_VPMPop));
    MinimalEventDistance = 40; % ms
    [~, EventIdx] = findpeaks(PSTH_VPMPop,'minpeakheight',MinimalEventHeight,'minpeakdistance',round(MinimalEventDistance/PSTHBinSize),'sortstr','descend');
    EventTimes = Edges_VPMPop(EventIdx); 
    SortedEventTimes = sort(EventTimes);
    
    FieldNames = cellfun(@(x) split(x,'_'), fieldnames(Model.model.fixed_variables),'uniformoutput',false);
    %% Analyse cortical (CTX) activity
    for pop = 1:length(CTXNames)
        PopNameLong = CTXNamesLong{pop};
        PopNameShort = CTXNames{pop};
        NeuronNumber = size(Model.(PopNameLong),2);
        
        %% Connectivity matrices
        NetCons = [];
        for fnms = 1:length(FieldNames)
            if strcmp(FieldNames{fnms, 1}{1, 1}, PopNameShort) && strcmp(FieldNames{fnms, 1}{end, 1}, 'netcon')
                NetCons.(cell2mat(join(FieldNames{fnms, 1},'_'))) = Model.model.fixed_variables.(cell2mat(join(FieldNames{fnms, 1},'_')));
            end
        end
        
        %% Average number of spikes per second, for every neuron
        SpikeTimes = cellfun(@(x) x(x>CutStartTime), Model.([PopNameLong '_spike_times']), 'uniformoutput', false);
        NumOfSpikes_PerSec = cellfun(@(x) numel(x)/SimulationTime_InSec, SpikeTimes, 'uniformoutput', false);
        
        %% Average ISI (inter-spike interval) per neuron
        AverageISI_InSec = cellfun(@(x) mean(diff((x)))/1000, SpikeTimes, 'uniformoutput', false);
        % if you want to get instantaneous FR, you have to calculate 1/AverageISI_InSec per neuron 
        
        %% Membrane potential (Vm) distribution
        FilteredVm = Model.(PopNameLong)((CutStartTime/SamplingPeriod)+1:end,:);
        VmDistributionBinEdges = [-100:0.5:50];
        VmDistribution = arrayfun(@(x) histcounts(FilteredVm(:,x), VmDistributionBinEdges), 1:NeuronNumber, 'uniformoutput', false)';
        % plot(VmDistributionBinEdges, [0 mean(cell2mat(VmDistribution))])
        
        %% full simulation cross correlation
        % single CTX cell PSTH with VPM population PSTH
            % CTX neuron activity and xcorr
        SingleNeuronPSTH_XCorr = [];
        SingleNeuronPSTH_Lags = [];
        for neuron = 1:NeuronNumber
            PSTH_SingleNeuron = plotPsthAndRaster('psth',SpikeTimes(neuron),PSTHBinSize,SimulationTime_InMs,0);
            [SingleNeuronPSTH_XCorr(neuron,:), SingleNeuronPSTH_Lags(neuron,:)] = xcorr(PSTH_VPMPop, PSTH_SingleNeuron);
        end
        
        % CTX population PSTH with VPM population PSTH
        PSTH_CTXPop = plotPsthAndRaster('psth',SpikeTimes,PSTHBinSize,SimulationTime_InMs,0);
        [CTX_PSTH_XCorr, CTX_PSTH_Lags] = xcorr(PSTH_VPMPop, PSTH_CTXPop);

        %% Event-based activity:
        % Number of VPM events
        EventNum = length(EventTimes);
        
        % Event-duration / Start+End times
        PreEvent = 20; % ms pre peak
        PostEvent = 60; % ms post peak
        PreBinNum = round(PreEvent/PSTHBinSize);
        PostBinNum = round(PostEvent/PSTHBinSize);
        
        VPM_EventData = [];
        for VPMEvent = 1:EventNum
            % Event start
            if (EventTimes(VPMEvent)-PreEvent) > 0
                VPM_EventData{VPMEvent}.StartTime = Edges_VPMPop(EventIdx(VPMEvent) - PreBinNum);
            else
                VPM_EventData{VPMEvent}.StartTime = 0;
                LeftPadding = zeros(1,abs(diff([EventIdx(VPMEvent) - PreBinNum, Edges_VPMPop(1)])));
            end
            % Event peak
            VPM_EventData{VPMEvent}.PeakTime = EventTimes(VPMEvent);
            
           % Event end
            if (EventTimes(VPMEvent)+PostEvent) < Edges_VPMPop(end)
                VPM_EventData{VPMEvent}.EndTime = Edges_VPMPop(EventIdx(VPMEvent) + PostBinNum);
            else
                VPM_EventData{VPMEvent}.EndTime = Edges_VPMPop(end);
                RightPadding = zeros(1,abs(diff([EventIdx(VPMEvent) + PostBinNum, Edges_VPMPop(end)])));
            end
            
            % Event raw activity (spike-times): from here I can make PSTH and then xcorr
            VPM_EventData{VPMEvent}.VPMSpikes = cellfun(@(x) x(x>VPM_EventData{VPMEvent}.StartTime & x<VPM_EventData{VPMEvent}.EndTime), FilteredVPM_Spikes, 'uniformoutput', false);
            VPM_EventData{VPMEvent}.CTXTime = cellfun(@(x) x(x>VPM_EventData{VPMEvent}.StartTime & x<VPM_EventData{VPMEvent}.EndTime), SpikeTimes, 'uniformoutput', false);
            
            % Number of neurons participating per event
            VPM_EventData{VPMEvent}.ParticipatingCTXNeurons = sum(cell2mat(cellfun(@(x) ~isempty(x), VPM_EventData{VPMEvent}.CTXTime, 'uniformoutput', false)));

            % Average ISI per neuron per event
            VPM_EventData{VPMEvent}.SpikeCount = cellfun(@(x) numel(x), VPM_EventData{VPMEvent}.CTXTime, 'uniformoutput', false);
            VPM_EventData{VPMEvent}.MeanISI = cellfun(@(x) mean(diff((x))), VPM_EventData{VPMEvent}.CTXTime, 'uniformoutput', false);

            % Spike latencies:
                % if the 1st the spike occurs within event_start+next_event_start 
                % and if this time is not within previous_event_start+previous_event_end, 
                % then I accept it
                % ^ the second condition should be only for L2/3
                
                % if the 2nd spike occurs within event_start+event_end 
                % and if this time is not within previous_event_start+previous_event_end
            VPM_EventData{VPMEvent}.FirstSpikeLatency = [];
            VPM_EventData{VPMEvent}.SecondSpikeLatency = [];
            
            SortedEventIdx = find(SortedEventTimes==EventTimes(VPMEvent));
            if SortedEventIdx < length(SortedEventTimes)
                NextEventPeak = SortedEventTimes(SortedEventIdx+1);
                TailToTheNextPeak = diff([VPM_EventData{VPMEvent}.PeakTime NextEventPeak]); % in ms
            else
                TailToTheNextPeak = PostEvent; % in ms
            end
            
            if  SortedEventIdx > 1
                PreviousEventPeak = SortedEventTimes(SortedEventIdx-1);
%                 TailFromPreviousToThisPeak = diff([PreviousEventPeak VPM_EventData{VPMEvent}.PeakTime]); % in ms
            else
                PreviousEventPeak = CutStartTime;
%                 TailFromPreviousToThisPeak = diff([VPM_EventData{VPMEvent}.PeakTime PreviousEventPeak]); % in ms
            end
            
            LeftSpikeBound = VPM_EventData{VPMEvent}.PeakTime; % in ms
            RightSpikeBound = VPM_EventData{VPMEvent}.PeakTime + TailToTheNextPeak; % in ms
            
            for neuron = 1:NeuronNumber
                OneNeuronSpikes = cell2mat(VPM_EventData{VPMEvent}.CTXTime(neuron));

                % in case there is more than 1 spike per
                % neuron, there is a chance that the first
                % spike occured before the event so its latency
                % is not registered because it does not meet
                % the criteria.. then, the second spike is
                % actually the first one
                OneNeuronSpikes = OneNeuronSpikes(OneNeuronSpikes>=LeftSpikeBound);
                
                if ~isempty(OneNeuronSpikes)
                    % First-spike latency
                    if length(OneNeuronSpikes) >= 1
                        FirstSpike = OneNeuronSpikes(1);
                        if FirstSpike >= LeftSpikeBound && FirstSpike < (RightSpikeBound - PreEvent)
                            FirstSpikeLatency = FirstSpike - VPM_EventData{VPMEvent}.PeakTime;
                            
                            % additional condition for spikes in L2/3
                            if any([contains(PopNameShort,'2'), contains(PopNameShort,'3')]) 
                                if  FirstSpike > (PreviousEventPeak + PostEvent)
                                    FirstSpikeLatency = FirstSpike - VPM_EventData{VPMEvent}.PeakTime;
                                else
                                    FirstSpikeLatency = [];
                                end
                            end
                        else
                            FirstSpikeLatency = [];
                        end

                        if length(OneNeuronSpikes) == 2
                            SecondSpike = OneNeuronSpikes(2);
                            % same conditions as for the first spike
                            if SecondSpike > LeftSpikeBound && SecondSpike < RightSpikeBound
                                SecondSpikeLatency = SecondSpike - VPM_EventData{VPMEvent}.PeakTime;
                                
                                % additional condition for spikes in L2/3
                                if any([contains(PopNameShort,'2'), contains(PopNameShort,'3')]) 
                                    if SecondSpike > (PreviousEventPeak + PostEvent)
                                        SecondSpikeLatency = SecondSpike - VPM_EventData{VPMEvent}.PeakTime;
                                    else
                                        SecondSpikeLatency = [];
                                    end
                                end
                            else
                                SecondSpikeLatency = [];
                            end
                        else % no second spike detected
                            SecondSpikeLatency = []; 
                        end
                    end
                else
                    FirstSpikeLatency = [];
                    SecondSpikeLatency = [];
                end
                VPM_EventData{VPMEvent}.FirstSpikeLatency{neuron} = FirstSpikeLatency;
                VPM_EventData{VPMEvent}.SecondSpikeLatency{neuron} = SecondSpikeLatency;
            end
        end
        
        % sample FR from nonETAs
        EventWindows = sort([0 EventTimes (EventTimes + PostEvent) Model.time(end)]);
        
        WindowPairs = reshape(EventWindows, [2, length(EventWindows)/2])';
        OutOfEvent_Start = WindowPairs(:,1);
        OutOfEvent_End = WindowPairs(:,2);
        % so, here we have an array of start-times (which are actually
        % end-times from previous events) coupled with end-times (which are
        % start-times from following events)
        
        % sample FR from ETAs
        EventWindows = sort([EventTimes (EventTimes + PostEvent)]);
        WindowPairs = reshape(EventWindows, [2, length(EventWindows)/2])';
        WithinEvent_Start = WindowPairs(:,1);
        WithinEvent_End = WindowPairs(:,2);
       
        DoThis = 1;
        if DoThis
            for nETA = 1:length(OutOfEvent_Start)
                nETA_activity = cellfun(@(x) x(x>OutOfEvent_Start(nETA) & x<OutOfEvent_End(nETA)), SpikeTimes, 'uniformoutput', false);
                if nETA<length(OutOfEvent_Start)
                    ETA_activity = cellfun(@(x) x(x>WithinEvent_Start(nETA) & x<WithinEvent_End(nETA)), SpikeTimes, 'uniformoutput', false);
                end
            end

            % randomly sampled spike times
            nETA_samples = length(cell2mat(nETA_activity'));
            ETA_samples = length(cell2mat(ETA_activity'));
            NumSamples = round((nETA_samples + ETA_samples) / 2);

            AllSpikeTimes = cell2mat(SpikeTimes');
            if NumSamples >= 1
                RTA_activity = AllSpikeTimes(randi(length(AllSpikeTimes), NumSamples,1));
            else
                RTA_activity = [];
            end
        else
            nETA_activity = [];
            ETA_activity = [];
            RTA_activity = [];
        end
        
        AnalysedResults = [table(...
            {lbl}, ... % Simulation label
            {SeedNum}, ... 
            {ScalingFactor}, ...
            {ScaledParameter}, ...
            {PopNameShort}, ... % Population name
            {SimulationTime_InMs}, ... % Simulation time in ms
            {SpikeTimes}, ... % Filtered spike times
            {NumOfSpikes_PerSec}, ... % Spikes/sec
            {AverageISI_InSec}, ... % Mean ISI [s]
            {VmDistribution}, ... % Membrane potential distribution
            {VmDistributionBinEdges}, ... % Membrane potential distribution edges
            {SingleNeuronPSTH_XCorr}, ... % XCorr between single-neuron PSTH and VPM PSTH
            {SingleNeuronPSTH_Lags}, ... % Lags for the above
            {CTX_PSTH_XCorr}, ... % XCorr between CTX population PSTH and VPM PSTH
            {CTX_PSTH_Lags}, ... % Lags for the above
            {EventNum}, ... % Number of VPM events
            {PSTHBinSize}, ... % Bin size for all PSTH calculations
            {VPM_EventData}, ... % Quantification of CTX activity during each event
            {nETA_activity},...
            {ETA_activity},...
            {RTA_activity},...
            {NetCons},...
            'VariableNames',{'SimulationLabel'; 'Seed'; 'ScalingFactor'; 'ScaledParameter'; 'Population'; 'SimulationDuration'; 'SpikeTimes'; 'SpikesPerSec'; 'MeanISI'; 'VmDistribution'; 'VmDistEdges';...
            'SingleNeuronPSTH_XCorr'; 'SingleNeuronPSTH_XCorrLags'; 'PopulationPSTH_XCorr'; 'PopulationPSTH_XCorrLags'; 'VPMEventNum'; 'PSTHBinSize'; 'VPMEventData';...
            'nETA_spikes'; 'ETA_spikes'; 'RTA_spikes'; 'ConnMx'}); AnalysedResults];    
    end
end

%% FUNCTIONS
function [out1, out2] = plotPsthAndRaster(WhichPlot, SpikeTimes, BinSize, EndTime, PlotOrNot, Col, timevec)
% This function plots psth and/or rasterplots with input arguments:
    % WhichPlot: 'psth', 'raster'
    % SpikeTimes are in miliseconds
    % BinSize of the histogram bins
    % endTime of the simulation (or x axis right limit)
    % plotOrNot is a binary input [0 or 1] depending on whether you just need data or also a plot output
    
% Optional
    % Col is a color to be used while plotting
    
    if nargin < 5
        disp('Not enough input arguments')
        return
    end
    
    if ~exist('Col','var')
       Col = 'k';
    end
    
    if iscell(SpikeTimes)
        if size(SpikeTimes,1) == length(SpikeTimes)
            SpikeTimes = SpikeTimes';
        end
        for c_idx = 1:length(SpikeTimes)
          this_cell = SpikeTimes{1,c_idx};
          if iscell(this_cell)
            this_cell = cell2mat(SpikeTimes{1,c_idx});
          end
          [psth_count(c_idx,:), psth_edges(c_idx,:)] = histcounts(this_cell,0:BinSize:EndTime);
        end
    else 
        fprintf('Only cell arrays with spike times per neuron are accepted for now. \n')
        return
    end

   switch WhichPlot
      case 'psth'
          out1 = [0 sum(psth_count,1) 0];
          out2 = [0 psth_edges(1,:)];
          if PlotOrNot
              if exist('timevec','var')
                  stairs(timevec,[sum(psth_count,1) 0], 'color',Col);
              else
                  stairs(psth_edges(1,:),[sum(psth_count,1) 0], 'color',Col);
              end
          end
      case 'raster'
          if PlotOrNot
              for c_idx = 1:length(SpikeTimes)
                  hold on;
                  this_cell = SpikeTimes{1,c_idx};
                  if iscell(this_cell)
                    this_cell = cell2mat(SpikeTimes{1,c_idx});
                  end
                  yVar = zeros([1, length(this_cell)]);
                  xVar = this_cell;
                  if size(xVar,1) ~= 1
                      xVar = xVar';
                  end
                  if size(yVar,1) ~= 1
                      yVar = yVar';
                  end
                  line([xVar; xVar], [yVar + (c_idx - 1); yVar + c_idx],'color',Col, 'linestyle','-','marker','none','linewidth',1)
              end
          end
          out1 = [];
          out2 = [];
    end 
end