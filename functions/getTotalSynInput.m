function tsyni = getTotalSynInput(model,popname,cell_idx)
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
    
    cutStartTime =  100/0.01 + 1; 
    for ii=1:length(presynPops)
        if contains(presynPops(ii),'I') 
            currentType = 'iGABAa';
            currentMonitor = 'IGABAa';
        else
            currentType = 'iAMPA';
            currentMonitor = 'IAMPA';
        end
        if contains(popname,'L4') && ii==3
            synCurrent_preCut = sum(model.([char(popname) '_' char(presynPops(ii)) '_' currentType '_' currentMonitor])(1:cutStartTime-1,cell_idx)) + ...
                sum(model.([char(popname) '_' char(presynPops(ii+1)) '_' currentType '_' currentMonitor])(1:cutStartTime-1,cell_idx));
            synCurrent_postCut = sum(model.([char(popname) '_' char(presynPops(ii)) '_' currentType '_' currentMonitor])(cutStartTime:end,cell_idx)) + ...
                sum(model.([char(popname) '_' char(presynPops(ii+1)) '_' currentType '_' currentMonitor])(cutStartTime:end,cell_idx));
            
            tsyni.VPM.settlingTime = synCurrent_preCut;
            tsyni.VPM.simulationTime = synCurrent_postCut;
        elseif contains(popname,'L4') && ii==4
            ...
        else
            synCurrent_preCut = sum(model.([char(popname) '_' char(presynPops(ii)) '_' currentType '_' currentMonitor])(1:cutStartTime-1,cell_idx));
            synCurrent_postCut = sum(model.([char(popname) '_' char(presynPops(ii)) '_' currentType '_' currentMonitor])(cutStartTime:end,cell_idx));
                    
            tsyni.(presynPops(ii)).settlingTime = synCurrent_preCut;
            tsyni.(presynPops(ii)).simulationTime = synCurrent_postCut;
        end
    end
end