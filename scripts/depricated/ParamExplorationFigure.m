addpath(genpath('/Users/teatompos/Desktop/Matlab_scripts'))
addpath(genpath('/Users/teatompos/Library/CloudStorage/GoogleDrive-t.tompos@neurophysiology.nl/.shortcut-targets-by-id/1f3yqxuDDFP_-nrGlSKbqcVC3RJJZZWTiWyK5Kl5x65maOi21i1OEIloSOR_UsCfsdiH2FbZf/Lab notebook Tea Tompoš/Thesis/Chapter 2 (The model)/Simulations/ParamSensitivity'))

clear all; close all; 

SaveResults = 1;
SaveFigures = 1;

synScalingFactors = [0.4:0.2:2.5];
pconnScalingFactors = [0.4:0.2:2.5];

seedNum = 1:5;
repNum = 5;

%% Check if files were simulated
cd '/Users/teatompos/Library/CloudStorage/GoogleDrive-t.tompos@neurophysiology.nl/.shortcut-targets-by-id/1f3yqxuDDFP_-nrGlSKbqcVC3RJJZZWTiWyK5Kl5x65maOi21i1OEIloSOR_UsCfsdiH2FbZf/Lab notebook Tea Tompoš/Thesis/Chapter 2 (The model)/Simulations/ParamSensitivity'
fileNames = dir;
fileNamesExist = [];
if length(seedNum)>1
    for s_num = 1:length(seedNum)
        fileIdx = cell2mat(cellfun(@(x) endsWith(x,sprintf('_seed%i.mat', s_num)), {fileNames.name}, 'uniformoutput', false))';
        fileNamesExist = [{fileNames(fileIdx).name}'; fileNamesExist];
    end
else
    fileIdx = cell2mat(cellfun(@(x) endsWith(x,sprintf('_seed%i.mat', seedNum)), {fileNames.name}, 'uniformoutput', false))';
    fileNamesExist = {fileNames(fileIdx).name}';
end

% files that have to be generated
fileNamesAll = [];
for s_num = 1:length(seedNum)
    progress = [];
    missing = [];
    s = seedNum(s_num);
    for synStr = 1:length(synScalingFactors)
        for pconn = 1:length(pconnScalingFactors)
            for r = 1:repNum
                NetSaveName = sprintf('N4_syn%i_pconn%i_rep%i_seed%i', round(synScalingFactors(synStr)*100), round(pconnScalingFactors(pconn)*100), r, s);
                fileNamesAll{end+1} = NetSaveName;
                
                if any(contains(fileNamesExist, NetSaveName))
                    progress(end+1) = 1;
                else
                    progress(end+1) = 0;
                    missing{end+1} = NetSaveName;
                end
            end
        end
    end
    fprintf('%0.2f%% of simulations are complete for seed %i.\n', (length(find(progress>0))/length(progress))*100,s)
    if ~isempty(missing) % show which files are missing
        missing';
    end
end

%% Check if simulations are analysed, if not -- then analyse the missing ones
seedNum = 3;
repNum = 1;
for s_num = 1:length(seedNum)
    s = seedNum(s_num);
    for synStr = 1:length(synScalingFactors)
        fails = 0;
        analysisFileName1 = sprintf('SNN_analysed_synStr%i_seed%i', round(synScalingFactors(synStr)*100), s);
        analysisFileName2 = sprintf('SNN_analysed_synStr%i', round(synScalingFactors(synStr)*100));
        try
            load(analysisFileName1)
        catch ME
            fprintf('Error for %s: %s.\n', analysisFileName1, ME.message)
            fails = fails + 1;
        end
        try
            SNN_analysed = [SNN_analysed; load(analysisFileName2)];
        catch ME
            fprintf('Error for %s: %s.\n', analysisFileName2, ME.message)
            fails = fails + 1;
        end
        
        if fails < 2
            analysedSims = SNN_analysed.label;
            fprintf('There are some analysed files... \n')
            analysedSims'
        elseif fails == 2
            analysedSims = 'none';
            fprintf('There are no analysed files... \n')
        end
            
        for pconn = 1:length(pconnScalingFactors)
            simLabel = sprintf('synStr_%i_pconn_%i_seed%i',round(synScalingFactors(synStr)*100), round(pconnScalingFactors(pconn)*100), s);
            if ~any(contains(analysedSims, simLabel))
                for r = 1:repNum
                    NetSaveName = sprintf('N4_syn%i_pconn%i_rep%i_seed%i', round(synScalingFactors(synStr)*100), round(pconnScalingFactors(pconn)*100), r, s);
                    load(NetSaveName)
                    fprintf('Analysing %s... \n', NetSaveName)
                    SNN_analysed = [SNN_analysed; analyzeSim(SpikingNetwork, [simLabel 'rep' num2str(r)], 0, 'threshold')];
                end
            end
        end
        fprintf('All synStr%s analysed (seed %i).\n', strrep(num2str(synStr),'.','', s))
        
        newAnalysedSim = SNN_analysed.label;
        if length(newAnalysedSim) > length(analysedSim)
            fprintf('Saving the new analysis..\n')
            save(analysisFileName1, 'SNN_analysed', '-v7.3')
        else
            fprintf('Nothing new was added to SNN_analysed so I will not save it.\n')
        end
    end
end

%% Further analysis
s_num = 4;

warning('off')
CortPopColorTags = [0.0000 0.4470 0.7410;... % blue,   l2E
    0.8500 0.3250 0.0980;... % orange,    l2I
    0.9290 0.6940 0.1250;... % yellow,    l3E
    0.4940 0.1840 0.5560;... % purple,    l3I
    0.4660 0.6740 0.1880;... % green,     l4E
    0.3010 0.7450 0.9330;... % light blue,l4I
    0.6350 0.0780 0.1840]; 

pops = ['L2E';'L2I';'L3E';'L3I';'L4E';'L4I'];
quant = [];
fig_count = 1;
row_num = 4;
% fig_num = round(length(synScalingFactors)/row_num);

for synStr = synScalingFactors
%     if any(ismember(fig_count, [1:row_num:fig_num*row_num]))
%         figure; clf; tiledlayout(row_num,length(pconnScalingFactors))
%     end
%     fig_count = fig_count+1;
    load(sprintf('SNN_analysed_synStr%s', strrep(num2str(synStr),'.','')))
    
    binSize_events=10;
    a=SNN_analysed.ETAs{1,1}.event_times{1, 1};
    event_start = (diff(a([1 2]))/binSize_events)-3;
    x_vec = ([1:(length(SNN_analysed.ETAs{1,1}.vpm_events{1, 1}))]-event_start)*binSize_events;
    fr_converter1 = (binSize_events/1000); % turn bins into seconds
    for pconn = pconnScalingFactors
%         nexttile; hold on;
        for jj=1:length(pops)
            label = ['synStr_' strrep(num2str(synStr),'.','') '_pconn_' strrep(num2str(pconn),'.','') '_seed' num2str(s_num)];
            tbl = SNN_analysed(strcmp(SNN_analysed.pop, pops(jj,:)) & strcmp(SNN_analysed.label, label),:);
            
%             if strcmp(pops(jj,:),pops(1,:))
%                 fr_converter2 = length(SNN_analysed.spike_times{7,1});
%                 data = cell2mat(cellfun(@(x) cell2mat(x.vpm_events'), tbl.ETAs,'uniformoutput',false));
%                 stdShade((data./fr_converter1)/fr_converter2,0.18,CortPopColorTags(end,:),x_vec) % vpm event
%                 plot([0 0], ylim,':k','HandleVisibility',"off")
%             end
            
            if ~isempty(tbl)
                %%% PLOT
%                 nexttile; hold on;
                fr_converter2 = length(tbl.spike_times{1,1}); % average over number of neurons N
                data = cell2mat(cellfun(@(x) cell2mat(x.ctx_events'), tbl.ETAs,'uniformoutput',false));
%                 stdShade((data./fr_converter1)/fr_converter2,0.18,CortPopColorTags(jj,:),x_vec) % ctx event
                
%                 xlabel('time (ms)')
%                 ylabel('(spikes/s) / N')
%                 xax = gca;
%                 xax.XTickMode = 'manual';
%                 xax.XTick = x_vec;
%                 xax.XTickLabel = x_vec;
%                 title(sprintf('%s; %0.1f; %0.1f',pops(jj,:), synStr, pconn))
                
                %%% RATIOMETRIC QUANTIFICATION
                data_fr = (data./fr_converter1)/fr_converter2;
                
                dff_max = []; dff_sig = [];
                for dtfr = 1:size(data_fr,1)
                    sig = data_fr(dtfr,:);
                    f0 = sig(1);
                    f1 = max(sig);
                    dff_max(end+1) = (f1-f0)/(f1+f0); % avoiding zero-division
                    dff_sig(end+1,:) = (sig-f0)./(sig+f0);
                end
                
%                 nexttile; hold on;
%                 stdShade(dff_sig,0.18,CortPopColorTags(jj,:),x_vec)
%                 
%                 xlabel('time (ms)')
%                 ylabel('FR change (%)')
%                 title(sprintf('%s; %0.1f; %0.1f',pops(jj,:), synStr, pconn))
%                 xax = gca;
%                 xax.XTickMode = 'manual';
%                 xax.XTick = x_vec;
%                 xax.XTickLabel = x_vec;
                                
                %%% FUNCTION FITTING
                upsample_factor = 10;
                signal = resample(mean(data_fr,1), upsample_factor, 1); % upsampled mean
                time = 1:length(signal); % i can not use x_vec because it has negative values
                
                % Define the functions
                gamma_func = @(b,time) double(b(1).*gampdf(time,b(2),max(b(3),eps)));
                bimodal_func = @(b,time) double(b(1).*exp(-((time-b(2)).^2)./(2*b(3).^2)) + b(4).*exp(-((time-b(5)).^2)./(2*b(6).^2)));
                beta_func = @(b,time) double(b(1).*betapdf(time,b(2),b(3)));
                lognormal_func = @(b,time) double(b(1).*lognpdf(time,b(2),b(3)));
                gaussian_func = @(b,time) double(b(1).*exp(-((time-b(2)).^2)./(2*b(3).^2)));
                
                % Initial guess
                beta0 = double([max(signal), mean(time), std(time), max(signal), mean(time), std(time)]); %  six parameters: a1, b1, c1, a2, b2, and c2. Here’s what they represent:
                % a1 and a2: These are the amplitudes of the two Gaussians.
                % b1 and b2: These are the means of the two Gaussians.
                % c1 and c2: These are the standard deviations of the two Gaussians.
                
                % Fit the functions
                models = {gamma_func, bimodal_func, beta_func, lognormal_func, gaussian_func};
                model_names = {'Gamma', 'Bimodal', 'Beta', 'Lognormal','Gaussian'};
                for ftype = 1:length(models)
                    model = models{ftype};
                    func_type = model_names{ftype};
                    if ~strcmp(func_type,'Bimodal')
                        beta0 = beta0(1:3);
                    end
                    try
                        [beta,R,~,~,MSE] = nlinfit(time, signal, model, beta0);
                        
                        SSE = sum(R.^2); % lower SSE, better fit; SSE provides an absolute measure of the discrepancy between the data and the model
                        RMSE = sqrt(MSE); % lower RMSE, better fit; RMSE provides a relative measure of fit, in units of the variable being predicted
                        
                        %%% get amplitude (peak of the distribution)
                        amplitude = beta(1);
                        
                        %%% get latency (mean of the distribution)
                        latency = beta(2); % the location of the center of the peak of the distribution
                        
                        %%% get slope (stdev of the distribution)
                        slope = beta(3); % smaller sigma/stdev means steeper slope (because the distribution is narrow)
                        
                        quant = [quant; table({pops(jj,:)},{pconn},{synStr},{func_type},{SSE},{RMSE},{amplitude},{latency},{slope},{dff_max},{dff_sig},...
                            'variablenames',{'pop','pconn','synstr','fit','SSE','RMSE','amplitude','latency','slope','df/f_max','df/f_sig'})];
                    catch ME
                        disp([func_type ' fitting did not work in ' label ' for ' pops(jj,:)])
                    end
                end
            else
                disp(['Could not find: ' label ' for ' pops(jj,:)])
                
                SSE = NaN;
                RMSE = NaN;
                amplitude = NaN;
                latency = NaN;
                slope = NaN;
                dff_max = NaN;
                dff_sig = NaN;
                func_type = '';
            end
%             if synStr == synScalingFactors(1) && pconn == pconnScalingFactors(end) && strcmp(pops(jj,:),pops(end,:))
%                 ax=get(gca,'Children');
%                 legend((ax(1:2:end)),[flipud(pops);'VPM'])
%             end
            
            quant = [quant; table({pops(jj,:)},{pconn},{synStr},{func_type},{SSE},{RMSE},{amplitude},{latency},{slope},{dff_max},{dff_sig},...
                'variablenames',{'pop','pconn','synstr','fit','SSE','RMSE','amplitude','latency','slope','df/f_max','df/f_sig'})];
        end
%         xlabel('time (ms)')
%         ylabel('(spikes/s) / N')
%         title(sprintf('%0.1f; %0.1f', synStr, pconn))
    end
end

if SaveResults
    cd '/Users/teatompos/Library/CloudStorage/GoogleDrive-t.tompos@neurophysiology.nl/.shortcut-targets-by-id/1f3yqxuDDFP_-nrGlSKbqcVC3RJJZZWTiWyK5Kl5x65maOi21i1OEIloSOR_UsCfsdiH2FbZf/Lab notebook Tea Tompoš/Thesis/Chapter 2 (The model)/Analysis'
    save(['ParamSensitivity_seed' num2str(s_num)], 'quant','-v7.3')
end


%% PLOT QUANTIFICATION
cd '/Users/teatompos/Library/CloudStorage/GoogleDrive-t.tompos@neurophysiology.nl/.shortcut-targets-by-id/1f3yqxuDDFP_-nrGlSKbqcVC3RJJZZWTiWyK5Kl5x65maOi21i1OEIloSOR_UsCfsdiH2FbZf/Lab notebook Tea Tompoš/Thesis/Chapter 2 (The model)/Analysis'

clear all;
load(['ParamSensitivity_seed' num2str(s_num)])
pops = ['L2E';'L2I';'L3E';'L3I';'L4E';'L4I'];
synScalingFactors = [0.4:0.2:2.5];
pconnScalingFactors = [0.4:0.2:2.5];

figure(1); clf; tiledlayout(2,3)
figure(2); clf; tiledlayout(2,3)
figure(3); clf; tiledlayout(2,3)
figure(4); clf; tiledlayout(2,3)

for jj = [1 3 5 2 4 6]
    M_amp = []; M_lat = []; M_slo = [];
    axis_max = [30 80 50 100;... % L2E
        30 80 50 100;... % L2I
        50 80 50 100;... % L3E
        50 80 50 100;... % L3I
        150 80 50 100;... % L4E
        150 80 50 100];  % L4I
    pop = pops(jj,:);
    for synStr = synScalingFactors
        for pconn = pconnScalingFactors
            tbl = quant(strcmp(quant.pop, pop) & cell2mat(quant.synstr)==synStr & cell2mat(quant.pconn)==pconn, :);
            if ~isempty(tbl.fit{1})
                %%% FIND THE FIT WITH THE SMALLEST ERROR AND PLOT IT
                best_fit = tbl(strcmp(tbl.pop, pops(jj,:)),:);
                best_fit_idx = find(cell2mat(best_fit.RMSE)==min(cell2mat(best_fit.RMSE)),1,'first');
                
                M_amp(find(synScalingFactors==synStr),find(pconnScalingFactors==pconn)) = best_fit.amplitude{best_fit_idx};    % matrix [synstr x pconn]
                M_lat(find(synScalingFactors==synStr),find(pconnScalingFactors==pconn)) = best_fit.latency{best_fit_idx};
                M_slo(find(synScalingFactors==synStr),find(pconnScalingFactors==pconn)) = best_fit.slope{best_fit_idx};
                
                M_f1f0(find(synScalingFactors==synStr),find(pconnScalingFactors==pconn)) = nanmean(best_fit.('df/f_max'){best_fit_idx});
            else
                M_amp(find(synScalingFactors==synStr),find(pconnScalingFactors==pconn)) = NaN;    % matrix [synstr x pconn]
                M_lat(find(synScalingFactors==synStr),find(pconnScalingFactors==pconn)) = NaN;
                M_slo(find(synScalingFactors==synStr),find(pconnScalingFactors==pconn)) = NaN;  
                M_f1f0(find(synScalingFactors==synStr),find(pconnScalingFactors==pconn)) = NaN;
            end
        end
    end
    %%% AMPLITUDE
    figure(1);nexttile;
    h = pcolor(smoothdata(M_amp, 'movmedian', 3.5));
    set(h, 'EdgeColor', 'none', 'AlphaData', ~isnan(M_amp), 'facecolor', 'interp');
    set(gca,'ydir','reverse')
    cb = colorbar;
    colormap(turbo(256));
    caxis([0,axis_max(jj,1)]);
    title(pop)
    ylabel(cb,'Amplitude (Spikes/s / N)')
    
    xax = gca;
    xax.XTickMode = 'manual';
    xax.XTick = 1:length(pconnScalingFactors);
    xax.XTickLabel = pconnScalingFactors; 
    
    xax.YTickMode = 'manual';
    xax.YTick = 1:length(synScalingFactors);
    xax.YTickLabel = synScalingFactors; 
    xlabel('PConn scale')
    ylabel('SynStr scale')
    
    %%% LATENCY
    figure(2);nexttile;
    h = pcolor(smoothdata(M_lat, 'movmedian', 3.5));
    set(h, 'EdgeColor', 'none', 'AlphaData', ~isnan(M_amp), 'facecolor', 'interp');
    set(gca,'ydir','reverse')
    cb = colorbar;
    colormap(turbo(256));
    caxis([0,axis_max(jj,2)]);
    title([pop '; Latency (ms)'])
    ylabel(cb, 'Peak response latency (ms)')
    
    xax = gca;
    xax.XTickMode = 'manual';
    xax.XTick = 1:length(pconnScalingFactors);
    xax.XTickLabel = pconnScalingFactors; 
    
    xax.YTickMode = 'manual';
    xax.YTick = 1:length(synScalingFactors);
    xax.YTickLabel = synScalingFactors; 
    xlabel('PConn scale')
    ylabel('SynStr scale')
    
    %%% SLOPES
    figure(3);nexttile;
    h = pcolor(smoothdata(M_slo, 'movmedian', 3.5));
    set(h, 'EdgeColor', 'none', 'AlphaData', ~isnan(M_amp), 'facecolor', 'interp');
    set(gca,'ydir','reverse')
    colormap(turbo(256));
    cb = colorbar;
    caxis([0,axis_max(jj,3)]);
    title(pop)
    ylabel(cb,'Slope')
    
    xax = gca;
    xax.XTickMode = 'manual';
    xax.XTick = 1:length(pconnScalingFactors);
    xax.XTickLabel = pconnScalingFactors; 
    
    xax.YTickMode = 'manual';
    xax.YTick = 1:length(synScalingFactors);
    xax.YTickLabel = synScalingFactors; 
    xlabel('PConn scale')
    ylabel('SynStr scale')
    
    %%% f1/f0
    figure(4);nexttile;
    h = pcolor(smoothdata(M_f1f0*100, 'movmedian', 3.5));
    set(h, 'EdgeColor', 'none', 'AlphaData', ~isnan(M_f1f0), 'facecolor', 'interp');
    set(gca,'ydir','reverse')
    colormap(turbo(256));
    cb = colorbar;
    caxis([0,axis_max(jj,4)]);
    title(pop)
    ylabel(cb,'FR change (%)')
    
    xax = gca;
    xax.XTickMode = 'manual';
    xax.XTick = 1:length(pconnScalingFactors);
    xax.XTickLabel = pconnScalingFactors; 
    
    xax.YTickMode = 'manual';
    xax.YTick = 1:length(synScalingFactors);
    xax.YTickLabel = synScalingFactors; 
    xlabel('PConn scale')
    ylabel('SynStr scale')
end

if SaveFigures
    cd '/Users/teatompos/Library/CloudStorage/GoogleDrive-t.tompos@neurophysiology.nl/.shortcut-targets-by-id/1f3yqxuDDFP_-nrGlSKbqcVC3RJJZZWTiWyK5Kl5x65maOi21i1OEIloSOR_UsCfsdiH2FbZf/Lab notebook Tea Tompoš/Thesis/Chapter 2 (The model)/Figures'
    fig_labels = {'amp'; 'lat'; 'slope'; 'ratio'};
    for fg =1:4
        figure(fg);
        f=gcf;
        set(f, 'paperorientation', 'landscape', 'renderer','painter')
        print(f,['AMPAx2_heatmap_' fig_labels{fg}], '-dpdf','-bestfit')
        close(f)
    end
end