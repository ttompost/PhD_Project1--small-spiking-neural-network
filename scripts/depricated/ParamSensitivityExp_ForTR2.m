cd '/Users/teatompos/Library/CloudStorage/GoogleDrive-t.tompos@neurophysiology.nl/.shortcut-targets-by-id/1f3yqxuDDFP_-nrGlSKbqcVC3RJJZZWTiWyK5Kl5x65maOi21i1OEIloSOR_UsCfsdiH2FbZf/Lab notebook Tea Tompoš/Thesis/Chapter 2 (The model)/Simulations/ParamSensitivity'
fileNames = dir;
fileIdx1 = cell2mat(cellfun(@(x) startsWith(x,'N4'), {fileNames.name}, 'uniformoutput', false))';
fileIdx2 = cell2mat(cellfun(@(x) startsWith(x,'SNN'), {fileNames.name}, 'uniformoutput', false))';
fileNamesExist = {fileNames(fileIdx1).name}';
analysisNamesExist = {fileNames(fileIdx2).name}';


synScalingFactors = [0.4:0.2:2.5];
pconnScalingFactors = [0.4:0.2:2.5];

seed_num = 1:5;
repetitions = 1:5;

Ecells_percentage = 0.85;
Icells_percentage = 0.15;


tstop = 600;

isLif = 0;
for s_num = seed_num
    for synStr = 4:length(synScalingFactors)
        SaveName = sprintf('SNN_analysed_synStr%s', strrep(num2str(synScalingFactors(synStr)),'.',''));
        if any(contains(analysisNamesExist, [SaveName '.mat']))
            fprintf('SynStr %0.2f is fully analysed.\n',synScalingFactors(synStr))
        else
            SNN_analysed = [];
            for pconn = 1:length(pconnScalingFactors)
                synParams = struct('synstr',[synScalingFactors(synStr) synScalingFactors(synStr)], ...
                    'pconn', [pconnScalingFactors(pconn) pconnScalingFactors(pconn)]);
                
                label = ['synStr_' strrep(num2str(synScalingFactors(synStr)),'.','') '_pconn_' strrep(num2str(pconnScalingFactors(pconn)),'.','') '_seed' num2str(s_num)];
                for r = repetitions
                    NetSaveName = sprintf('N4_syn%s_pconn%s_rep%i_seed%i', num2str(synScalingFactors(synStr)*100), num2str(pconnScalingFactors(pconn)*100), r, s_num);
                    if ~any(contains(fileNamesExist, [NetSaveName '.mat']))
                        fprintf('Simulating %s.\n',NetSaveName)
                        SpikingNetwork = SimN4(Ecells_percentage, Icells_percentage, [], synParams, [], [], 'synthetic', s_num, tstop);
                    
%                         save(NetSaveName,'SpikingNetwork','-v7.3')
                        ParallelSave(NetSaveName,SpikingNetwork)
                    else
                        fprintf('Loading %s.\n',NetSaveName)
%                         load([NetSaveName '.mat'])

                    end
                    fprintf('Analysing %s.\n',NetSaveName)
                    SNN_analysed = [SNN_analysed; analyzeSim(SpikingNetwork, label, isLif, 'threshold')];
                    save(SaveName,'SNN_analysed','-v7.3')
                end
            end
            save(SaveName,'SNN_analysed','-v7.3')
%             ParallelSave(SaveName,SNN_analysed)
        end
    end
end