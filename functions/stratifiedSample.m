function [sampledData1] = stratifiedSample(data1, sampleSize)

% Input:
%   data1: The original data set with 9k elements (assumed to be a vector)
%   sampleSize: The desired size of the subsample (e.g., size of data2)

% Output:
%   sampledData1: A random subsample of data1 with size 'sampleSize' 
%                  preserving the distribution

    % Check input validity
    if sampleSize > length(data1)
      error('Sample size cannot be larger than the original data size');
    end

    % Count elements in each unique value group (stratum)
    [counts, edges1] = histcounts(data1, sampleSize);
    [binIdices, edges2] = discretize(data1, sampleSize);
    if ~sum(edges1 == edges2) == length(edges1)
        error('Histcounts() and Discretize() outputs do not match.');
    end
    
    % Calculate proportion of data points in each stratum
    proportions = counts / sum(counts);

    % Pre-allocate for sampled data
    sampledData1 = [];

    % Loop through each unique value (stratum)
    for idx = 1:sampleSize
      % Calculate the number of samples to take from this stratum
      if proportions(idx) > 0
          stratumSampleSize = round(sampleSize * proportions(idx));

          % Sample 'stratumSampleSize' elements from the current stratum (unique value)
          stratumIndices = find(binIdices == idx);
          sampledData = data1(stratumIndices(1:stratumSampleSize));
          
          sampledData1(end+1:end+length(sampledData),1) = sampledData;
      end
    end
    
    % rounding error produces results that can be way off the desired
    % sampleSize, so here i randomly sample values from data1 to fill that
    % gap
    while length(sampledData1) < sampleSize || length(sampledData1) > sampleSize
        if length(sampledData1) < sampleSize
            randomSamples = randi(sampleSize, sampleSize - length(sampledData1),1);
            sampledData1(end+1:end+length(randomSamples),1) = data1(randomSamples);
        elseif length(sampledData1) > sampleSize
            randomSamples = randi(length(sampledData1), length(sampledData1) - sampleSize,1);
            sampledData1(randomSamples) = [];
        end
    end
end