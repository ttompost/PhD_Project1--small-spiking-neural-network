function tsyni = getTotalSynInput(model,popname,cell_idx,timewindows)
    % tsyni = total synaptic current
    
    if contains(popname,'L2')
        presynPops = ["L2E", "L2I","L3E","L3I"];
    elseif contains(popname,'L3E')
        presynPops = ["L3E","L3I","L4E", "L4I"];
    elseif contains(popname,'L3I')
        presynPops = ["L3E","L3I","L4E"];
    elseif contains(popname,'L4')
        presynPops = ["L4E","L4I","VPMcorr", "VPMuncorr"];
    end
    
    % to get the total synaptic current per cell in a population, it is
    % enough to work with (e.g.) model.L2E_L3E_iAMPA_IAMPA matrix because
    % SpikingNetwork.L2E_L3E_iAMPA_IAMPA(:,cell_idx) contains AMPAergic
    % activity from all L3E presynaptic partners to the postsynaptic L2E
    % cell_idx cells
    
    dt = 0.01;
    tsyni = [];
    for ii = 1:length(presynPops)
        if contains(presynPops(ii),'I') % inhibitory presynapses
            currentType = 'iGABAa';
            currentMonitor = 'IGABAa';
        else % excitatory presynapses
            currentType = 'iAMPA';
            currentMonitor = 'IAMPA';
        end
        if exist('timewindows','var')
            for tw = 1:size(timewindows,1)
                winStart = find(round(model.time,2) == round(timewindows(tw, 1),2));
                winEnd = find(round(model.time,2) == round(timewindows(tw, 2),2));
                if contains(popname,'L4') && ii==3  % VPM presynapses for L4
                    tsyni.VPM.timeWindows{tw} = (model.([char(popname) '_' char(presynPops(ii)) '_' currentType '_' currentMonitor])(winStart:winEnd,cell_idx)) + ...
                        (model.([char(popname) '_' char(presynPops(ii+1)) '_' currentType '_' currentMonitor])(winStart:winEnd,cell_idx));
                elseif contains(popname,'L4') && ii==4
                    ... % VPM will be detected twice, but it was already analysed in the upper condition
                else
                    tsyni.(presynPops(ii)).timeWindows{tw} = (model.([char(popname) '_' char(presynPops(ii)) '_' currentType '_' currentMonitor])(winStart:winEnd,cell_idx));
                end
            end
            % add full simulation in the end
            winStart = 1;
            winEnd = length(model.time);
            if contains(popname,'L4') && ii==3  % VPM presynapses for L4
                tsyni.VPM.fullSimulation{1} = (model.([char(popname) '_' char(presynPops(ii)) '_' currentType '_' currentMonitor])(winStart:winEnd,cell_idx)) + ...
                    (model.([char(popname) '_' char(presynPops(ii+1)) '_' currentType '_' currentMonitor])(winStart:winEnd,cell_idx));
            elseif contains(popname,'L4') && ii==4
                ... % VPM will be detected twice, but it was already analysed in the upper condition
            else
                tsyni.(presynPops(ii)).fullSimulation{1} = (model.([char(popname) '_' char(presynPops(ii)) '_' currentType '_' currentMonitor])(winStart:winEnd,cell_idx));
            end
             
        else
            tw = 1;
            winStart = 1;
            winEnd = length(model.time);
            if contains(popname,'L4') && ii==3  % VPM presynapses for L4
                tsyni.VPM.timeWindows{tw} = (model.([char(popname) '_' char(presynPops(ii)) '_' currentType '_' currentMonitor])(winStart:winEnd,cell_idx)) + ...
                    (model.([char(popname) '_' char(presynPops(ii+1)) '_' currentType '_' currentMonitor])(winStart:winEnd,cell_idx));
            elseif contains(popname,'L4') && ii==4
                ... % VPM will be detected twice, but it was already analysed in the upper condition
            else
                tsyni.(presynPops(ii)).timeWindows{tw} = (model.([char(popname) '_' char(presynPops(ii)) '_' currentType '_' currentMonitor])(winStart:winEnd,cell_idx));
            end
            
        end
    end
end