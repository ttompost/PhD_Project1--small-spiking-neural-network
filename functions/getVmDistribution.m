function [d, y] = getVmDistribution(x)
    BinSize = 2; %in mV
    d=[]; y=[];
    for ii=1:size(x,2)
        sig = x(:,ii);
        sig=sig(sig<0);% cut APs above 0mv
        [d(ii,:), y(ii,:)]= histcounts(sig,-100:BinSize:20);
    end
    tp = numel(x); % total number of points
    d=[sum(d) 0]/tp; % make a fraction
%     d=d*100; % make percentage
    y=y(ii,:);
end