function kl = KL_Div_hist(p,q,x)
% computes KL divergence between 2 histograms

% how well does q (chosen arbitrary distribution) approximate p (true distribution)

    % first divide every column by the total count to convert counts to proportions, so you could treat the histograms as probability mass
    p = p./sum(p);
    q = q./sum(q);
    
    % now we take every bin and get p(bin) * log (p(bin)/q(bin))
    for ii = 1:length(x)-1
        y(ii) = p(ii) * log10(p(ii)/q(ii));
    end
    
    % remove nans
    y(isnan(y)) = [];
    % remove inf
    y(y==Inf) = [];
    y(y==-Inf) = [];
    
    % get KL divergence
    kl = sum(y);
end