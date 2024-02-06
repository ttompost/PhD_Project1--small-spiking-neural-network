addpath(genpath('/Users/teatompos/Desktop/Matlab_scripts/'))
clear all; close all; 

% set the network seed (this should be based on the exclusion criteria experiment)
netSeed = randperm(100,15);

% set the scaling factor array
ScalingFactors = [12:-0.5:1.5 1:-0.05:0.01];
ChosenSFs = [1:-0.1:0.1];

figure; clf; hold on;
plot(1:length(ScalingFactors), ScalingFactors, '-k','linewidth', 1.5)
plot(1:length(ScalingFactors), ScalingFactors, 'ok','linewidth', 1.5)
plot(find(ismember(round(ScalingFactors,2), round(ChosenSFs,2))), ChosenSFs, 'or','linewidth', 1.5)

plot([0 length(ScalingFactors)],[1 1], ':r', 'linewidth', 1.5)
xlim([0 length(ScalingFactors)])

text(3,0.5,'Sparsification', 'fontsize',16)
text(3,1.5,'Densification', 'fontsize',16)

xlabel('Scaling factor (number)')
ylabel('Scaling factor (value)')
set(gca,'fontsize', 14)

varyType = 'probs';

projections = [ "all Th to L4","all L4 to L3","all L3 to L2"... % translaminar global
        "ThE to L4E","ThE to L4I","l4E to L3E","L4E to L3I","L4I to L3E","L3E to L2E","L3E to L2I","L3I to L2E","L3I to L2I"... % translaminar local
        "all (intra)layers E to E","all (intra)layers E to I","all (intra)layers I to E","all (intra)layers I to I"... % intralaminar global
        "L4E intralayer", "L4E to L4I", "L4I to L4E", "L4I intralayer", "L3E intralayer", "L3E to L3I", "L3I to L3E", "L3I intralayer", "L2E intralayer", "L2E to L2I", "L2I to L2E", "L2I intralayer"... % intralaminar local
               ];

repetitions = 10;

Ecells_percentage = 0.85;
Icells_percentage = 0.15;

synParams = struct('synstr',[1 1], 'pconn', [1 1]);

tstop = 1100;

% % how many workers are available for parallel processing
% gcp_props = gcp;
% gcp_props.NumWorkers

% % how many workers do you wanna use?
% workNum = 9;

cd '/Users/teatompos/Library/CloudStorage/GoogleDrive-t.tompos@neurophysiology.nl/.shortcut-targets-by-id/1f3yqxuDDFP_-nrGlSKbqcVC3RJJZZWTiWyK5Kl5x65maOi21i1OEIloSOR_UsCfsdiH2FbZf/Lab notebook Tea Tompoš/Thesis/Chapter 2 (The model)/Simulations/Sparsification'
fileNames = dir;
fileIdx = cell2mat(cellfun(@(x) startsWith(x,'N4'), {fileNames.name}, 'uniformoutput', false))';
fileNamesExist = {fileNames(fileIdx).name}';
% fileNamesExist = {'none'};

% % files that have to be generated
% fileNamesAll = [];
% for sp = 1:length(projections)
%     progress = [];
%     sparseProjection = projections(sp);
%     for sf = find(ScalingFactors==1):2:length(ScalingFactors) % only sparsify
%         for s_Idx = 1:length(netSeed)
%             NetSaveName = sprintf('N4_sparse%i_%s_sfIdx_%i_seed%i', sp, regexprep(sparseProjection, ' ', '_'), sf, s_Idx);
%             fileNamesAll{end+1} = NetSaveName;
%             
%             progress(end+1) = any(contains(fileNamesExist, NetSaveName));
%         end
%     end
%     fprintf('(%i) %s: %0.2f%% of simulations are complete.\n', sp, regexprep(sparseProjection, ' ', '_'), (length(find(progress>0))/length(progress))*100)
% end

cd '/Users/teatompos/Library/CloudStorage/GoogleDrive-t.tompos@neurophysiology.nl/.shortcut-targets-by-id/1f3yqxuDDFP_-nrGlSKbqcVC3RJJZZWTiWyK5Kl5x65maOi21i1OEIloSOR_UsCfsdiH2FbZf/Lab notebook Tea Tompoš/Thesis/Chapter 2 (The model)/Simulations/Sparsification'

failedToSave = [];
failedToSim = [];
for sp = 5:length(projections)
    sparseProjection = projections(sp);
    
    warning('off')
    for sf = find(ScalingFactors==1):2:length(ScalingFactors) % only sparsify
        sparseFactor = ScalingFactors(sf);
        sparse_info = [];
        
        sparse_info{1} = sp;  % which projection is being sparsified
        sparse_info{2} = varyType;          % which connectivity parameter is affected
        sparse_info{3} = sparseFactor;      % what scaling factor is applied
        
        for s_Idx = 1:length(netSeed)
            NetSaveName = sprintf('N4_sparse%i_%s_sfIdx_%i_seed%i', sp, regexprep(sparseProjection, ' ', '_'), sf, s_Idx);
            if ~any(contains(fileNamesExist, [NetSaveName '.mat']))
                try
                    fprintf('Simulating %s.', NetSaveName)
                    SpikingNetwork = SimN4(Ecells_percentage, Icells_percentage, [], synParams, [], [], 'synthetic', netSeed(s_Idx), tstop, sparse_info);
                    SpikingNetwork.SparseExpInfo = sparse_info;
                    SpikingNetwork.SeedNum = netSeed(s_Idx);      % which seed is used
                catch ME
%                     failedToSim{sp} = NetSaveName;
                    fprintf('Failed to simulate: %s', NetSaveName)
                end
                
                try
                    save(NetSaveName,'SpikingNetwork','-v7.3')
%                     ParallelSave(NetSaveName,SpikingNetwork)
                catch ME
%                     failedToSave{sp} = NetSaveName;
                    fprintf('Failed to save: %s', NetSaveName)
                end
            end
        end
    end
end
