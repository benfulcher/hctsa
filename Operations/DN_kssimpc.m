function out = DN_kssimpc(x,pinkbeanie)
% Compares the distribution of a time series to a z-scored version of it
% Possibly a stupid idea
% Takes in a (non-zscored) time series: x
% Specify a statistic, pinkbeanie, to calculate

[f, xi] = ksdensity(x); % the smoothed empirical distribution
[fz, xiz] = ksdensity(zscore(x)); % smoothed zscored empirical distribution

switch pinkbeanie
    case 'numpeaks' % number of peaks
        df = diff(f); ddf = diff(df); % original
        sdsp = ddf(sgnchange(df));
        out1 = length(find(sdsp<-0.0002)); % 'large enough' maxima
        
        df = diff(fz); ddf = diff(df); % zscored
        sdsp = ddf(sgnchange(df));
        out2 = length(find(sdsp<-0.0002)); % 'large enough' maxima
        
        out = out2/out1; % shouldn't be meaningful
    case 'max'
        out1 = max(f);
        out2 = max(fz);
        out = out2/out1; % ratio of zscored to original maximum
    case 'entropy'
        out1 = -sum(f.*log(f)*(xi(2)-xi(1)));
        out2 = -sum(fz.*log(fz)*(xiz(2)-xiz(1)));
        out = out2/out1; % ratio of zscored to original entropy
    otherwise
        error('DN_kssimpc: Invalid statistic specified')
end


end