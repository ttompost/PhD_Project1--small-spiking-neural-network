function estimated_mean = EstimateMean_FromDistribution(counts, edges)
    midpoints = (edges(1:end-1) + edges(2:end)) / 2;
    total_count = sum(counts);
    weighted_sum = sum(midpoints .* counts);
    estimated_mean = weighted_sum / total_count;
end
