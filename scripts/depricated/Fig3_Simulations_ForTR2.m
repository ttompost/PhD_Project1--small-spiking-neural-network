% Path to the code/scripts
addpath(genpath('/Users/teatompos/Desktop/git_files'))
addpath(genpath('/Users/teatompos/Desktop/Matlab_scripts/RU_1st_internship_Neurophysiology_Dept/DynaSim'))

% Saving names prefixes

SaveSim_pre = 'Fig3_Net_';
SaveAn_pre = 'Fig3_Analysis_';
SaveFig_pre = 'Fig3_';

% Path to save/retrieve results
ResultsPath = '/Users/teatompos/Library/CloudStorage/GoogleDrive-t.tompos@neurophysiology.nl/.shortcut-targets-by-id/1f3yqxuDDFP_-nrGlSKbqcVC3RJJZZWTiWyK5Kl5x65maOi21i1OEIloSOR_UsCfsdiH2FbZf/Lab notebook Tea Tompoš/Thesis/Chapter 2 (The model)/Simulations/ParamSensitivity/Fixed_VPM_input';
addpath(genpath(ResultsPath))

Simulate = 1;
SaveResults = 1;
Analyse = 1;

cd(ResultsPath)

% Scaling factors
synScalingFactors = [0.4:0.2:2.5];
pconnScalingFactors = [0.4:0.2:2.5];

% E-I ratio
Ecells_percentage = 0.85;
Icells_percentage = 0.15;

% Simulation time
tstop = 600;

% LIF neurons?
isLif = 0;

% Which seeds to use?
RandomSeeds = [11, 48, 23, 39, 83];
seeds = [];

% Repetition per seed
repetitions = 10;

% Stimulate different networks with identical thalamic input
if Simulate
    for s_num = RandomSeeds(3)
        for synStr = 1:length(synScalingFactors)
            failedStuff = [];
            newSim = 0;
            newAn = 0;
            
            SNN_analysed = [];

            for pconn = 1:length(pconnScalingFactors)
                synParams = struct('synstr',[synScalingFactors(synStr) synScalingFactors(synStr)], ...
                        'pconn', [pconnScalingFactors(pconn) pconnScalingFactors(pconn)]);

               rng('shuffle') % for every synstr-pconn pair, create an array of random ctx seeds used in simulations below
               CTX_seeds = randi(100,repetitions,1); 
               
               for r_num = 1:repetitions
                    NetSaveName = sprintf('%s_S%i_P%i_Rep%i_Seed%i', SaveSim_pre, round(synScalingFactors(synStr)*100), round(pconnScalingFactors(pconn)*100), r_num, s_num);


                    seeds = struct('ctx_seed', CTX_seeds(r_num),... % randomize ctx connectivity for every repetition
                            'vpm_seed', s_num); % fix the vpm activity across repetitions

                    SpikingNetwork = SimN4(Ecells_percentage, Icells_percentage, [], synParams, [], [], 'synthetic', seeds, tstop);
                    SpikingNetwork.Seeds = seeds;
                    fprintf('Simulated %s. (CTX seed: %i; VPM seed: %i) \n', NetSaveName, seeds.ctx_seed, seeds.vpm_seed)

                    newSim = 1;
                    newAn = 1;

                    if SaveResults
                        if newSim
                            cd(ResultsPath)
                            save(NetSaveName, 'SpikingNetwork', '-v7.3')
                        end
                    end

                    if Analyse
                        if newAn
                            SNN_analysed = [SNN_analysed; analyzeSim(SpikingNetwork, NetSaveName, isLif, 'threshold')];

                            % inermediate save
                            cd(AnalysisPath)
                            AnSaveName = sprintf('%s_S%i_Seed%i_intermediate', SaveAn_pre, round(synScalingFactors(synStr)*100), s_num);
                            save(AnSaveName, 'SNN_analysed', '-v7.3')
                        end
                    end
                    newSim = 0;
                    newAn = 0;
               end
            end
            if SaveAnalysis
                cd(AnalysisPath)
                AnSaveName = sprintf('%s_S%i_Seed%i', SaveAn_pre, round(synScalingFactors(synStr)*100), s_num);
                save(AnSaveName, 'SNN_analysed', '-v7.3')
                
                fprintf('Saved the following analysis: %s.\n', AnSaveName)
            end
        end
    end
end