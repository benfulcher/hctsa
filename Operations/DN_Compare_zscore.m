% DN_Compare_zscore
% 
% Compares the distribution of a time series to a z-scored version of it
% 
% INPUTS:
% x, a (not z-scored) time series
% 
% Outputs are ratios of features between the original and z-scored time series,
% including the number of peaks, the maximum, and the distributional entropy.
% 

function out = DN_Compare_zscore(x)
% Specify a statistic, whatfeature, to compare

[f, xi] = ksdensity(x); % the smoothed empirical distribution
[fz, xiz] = ksdensity(BF_zscore(x)); % smoothed z-scored empirical distribution

% 1. numpeaks

df = diff(f);
ddf = diff(df); % original
sdsp = ddf(BF_sgnchange(df));
out1 = sum(sdsp < -0.0002); % 'large enough' maxima

df = diff(fz);
ddf = diff(df); % zscored
sdsp = ddf(BF_sgnchange(df));
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