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
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = DN_Compare_zscore(x)

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


end