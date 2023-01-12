function [out1, out2] = plotPsthAndRaster(WhichPlot, SpikeTimes, BinSize, EndTime, PlotOrNot, Col)
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
    
    if iscell(SpikeTimes)
      for c_idx = 1:length(SpikeTimes)
          this_cell = SpikeTimes{1,c_idx};
          [psth_count(c_idx,:), psth_edges(c_idx,:)] = histcounts(this_cell,0:BinSize:EndTime);
      end
    else 
        fprintf('Only cell arrays with spike times per neuron are accepted for now. \n')
        return
    end

   switch WhichPlot
      case 'psth'
          out1 = [sum(psth_count,1) 0];
          out2 = psth_edges(1,:);
          if PlotOrNot
              if exist('Col','var')
                  stairs(psth_edges(1,:),[sum(psth_count,1) 0], 'color',Col);
              else
                  stairs(psth_edges(1,:),[sum(psth_count,1) 0]);
              end
          end
      case 'raster'
          if PlotOrNot
              for c_idx = 1:length(SpikeTimes)
                  hold on;
                  yVar = zeros([1, length(SpikeTimes{1,c_idx})]);
                  xVar = SpikeTimes{1,c_idx};
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