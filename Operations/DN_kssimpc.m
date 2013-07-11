function out = DN_kssimpc(x,whatfeature)
% Compares the distribution of a time series to a z-scored version of it
% Possibly a stupid idea
% Takes in a (non-zscored) time series: x
% Specify a statistic, whatfeature, to compare

[f, xi] = ksdensity(x); % the smoothed empirical distribution
[fz, xiz] = ksdensity(benzscore(x)); % smoothed z-scored empirical distribution

% 1. numpeaks

df = diff(f);
ddf = diff(df); % original
sdsp = ddf(sgnchange(df));
out1 = sum(sdsp < -0.0002); % 'large enough' maxima

df = diff(fz);
ddf = diff(df); % zscored
sdsp = ddf(sgnchange(df));
out2 = sum(sdsp < -0.0002); % 'large enough' maxima

out.numpeaks = out2/out1; % shouldn't be meaningful

% 2. Max
out1 = max(f);
out2 = max(fz);
out.max = out2/out1; % ratio of zscored to original maximum

% 3. Entropy
out1 = -sum(f.*log(f)*(xi(2)-xi(1)));
out2 = -sum(fz.*log(fz)*(xiz(2)-xiz(1)));
out.entropy = out2/out1; % ratio of z-scored to original entropy


% switch whatfeature
%     case 'numpeaks' % number of peaks
%     case 'max'
%     case 'entropy'
%     otherwise
%         error('Invalid statistic specified: should be ''numpeaks'', ''max'', or ''entropy''')
% end


end