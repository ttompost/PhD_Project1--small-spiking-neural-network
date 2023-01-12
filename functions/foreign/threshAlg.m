function [peaks,avgFilter,stdFilter] = threshAlg(y,lag,threshold,influence)
%%% source: 
%%% https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/54507329#54507329


%%% calculates peaks/signals based on the principle of dispersion: 
% if a new datapoint is a given x number of standard deviations away 
% from some moving mean, the algorithm signals (also called z-score). 

% input params
    % lag - the lag of the moving window 
        % e.g. a lag of 5 will use the last 5 observations to smooth the data
        % The more stationary your data is, the more lags you should include 
        % (this should improve the robustness of the algorithm). If your data 
        % contains time-varying trends, you should consider how quickly you want 
        % the algorithm to adapt to these trends. I.e., if you put lag at 10, 
        % it takes 10 'periods' before the algorithm's treshold is adjusted to 
        % any systematic changes in the long-term average. 
        
    % threshold - the z-score at which the algorithm detects peaks
        % e.g. a threshold of 3.5 will signal if a datapoint is 3.5 standard deviations 
        % away from the moving mean
        
    % influence - the influence (between 0 and 1) of new peak on the mean and standard deviation
        % e.g. an influence of 0.5 gives signals half of the influence that normal datapoints have; 
        % Likewise, an influence of 0 ignores signals completely for recalculating the new threshold;
        % An influence of 0 is therefore the most robust option (but assumes stationarity); 
        % putting the influence option at 1 is least robust. For non-stationary data, 
        % the influence option should therefore be put somewhere between 0 and 1.


% Initialise signal results
peaks = zeros(length(y),1);
% Initialise filtered series
filteredY = y(1:lag+1);
% Initialise filters
avgFilter(lag+1,1) = mean(y(1:lag+1));
stdFilter(lag+1,1) = std(y(1:lag+1));
% Loop over all datapoints y(lag+2),...,y(t)
for i=lag+2:length(y)
    % If new value is a specified number of deviations away
    if abs(y(i)-avgFilter(i-1)) > threshold*stdFilter(i-1)
        if y(i) > avgFilter(i-1)
            % Positive signal
            peaks(i) = 1;
        else
            % Negative signal
            peaks(i) = -1;
        end
        % Make influence lower
        filteredY(i) = influence*y(i)+(1-influence)*filteredY(i-1);
    else
        % No signal
        peaks(i) = 0;
        filteredY(i) = y(i);
    end
    % Adjust the filters
    avgFilter(i) = mean(filteredY(i-lag:i));
    stdFilter(i) = std(filteredY(i-lag:i));
end
% Done, now return results
end