function out = WL_modwtvar(y, wname, level)
% WL_modwtvar   Multiscale variance decomposition via the maximal overlap DWT.
%
% Decomposes the time series into octave-scale bands using the maximal overlap
% discrete wavelet transform (MODWT) and computes summary statistics on how
% variance is distributed across scales. Unlike the standard (decimated) DWT used
% elsewhere in this codebase, the MODWT is shift-invariant and its variance
% estimator uses only the coefficients unaffected by the boundary at each level
% (Percival and Walden, 2000): the mean square of the unaffected coefficients of
% each detail level, and the variance about the mean of those of the scaling band.
%
% ---INPUTS:
% y, the input time series
% wname, the mother wavelet, e.g., 'db3', 'sym2' (see Wavelet Toolbox
%        Documentation; default: 'db3')
% level, the level of wavelet decomposition (can be set to 'max' for the maximum
%        level supported by the series length and wavelet,
%        floor(log2((N - 1)/(L - 1) + 1)) with L the filter length; default: 5). A
%        level that is too large is reduced to the maximum, so that every level, and
%        the scaling band, has at least one coefficient unaffected by the boundary.
%
% ---OUTPUTS:
% scalingFrac, the fraction of variance in the lowest-frequency (scaling/trend)
%        band: the part of the signal not resolved by any detail level
% domlevel, the level (1 = highest frequency, ..., level+1 = scaling/trend band)
%        that carries the most variance
% decaySlope, the slope of log2(variance) against level across the detail bands: a
%        wavelet-based scaling exponent, analogous to a Hurst estimate but using the
%        MODWT's unbiased, boundary-corrected variance rather than an ad hoc
%        regression on raw coefficients
%
% ---REFERENCES:
% D.B. Percival and A.T. Walden, "Wavelet Methods for Time Series Analysis",
% Cambridge University Press, Cambridge (2000). DOI: 10.1017/CBO9780511841040
% (the MODWT is treated on pp. 159-205 and the wavelet variance on pp. 295-339).

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
%% Check that a Wavelet Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('wavelet_toolbox');

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
y = y(:);
N = length(y); % length of the time series

if nargin < 2 || isempty(wname)
	wname = 'db3'; % Daubechies wavelet filter
end
if nargin < 3 || isempty(level)
	level = 5; % level of wavelet decomposition
end

% Deepest level whose scale has at least one boundary-unaffected coefficient, so that
% every retained level (and the scaling band) has an unbiased variance estimate:
% (2^j - 1)(L - 1) + 1 <= N, with L the wavelet filter length
filtLen = numel(wfilters(wname));
maxLevelAllowed = floor(log2((N - 1) / (filtLen - 1) + 1));
if strcmp(level, 'max')
	level = maxLevelAllowed;
end
if maxLevelAllowed < level
	fprintf(1, 'Chosen level (%u) is too large for a series of length N = %u\n', level, N);
	level = maxLevelAllowed;
	fprintf(1, 'Using a level of %u instead.\n', level);
end
if level < 1
	error('Time series too short for a MODWT decomposition');
end

% ------------------------------------------------------------------------------
%% MODWT and its unbiased multiscale variance decomposition
% ------------------------------------------------------------------------------
w = modwt(y, wname, level); % (level+1) x N: levels 1..level (detail), then scaling/trend

% Unbiased variance by level (Percival and Walden, 2000): level j has
% L_j = (2^j - 1)(L - 1) + 1 filter taps, so its last M_j = N - L_j + 1 coefficients
% are unaffected by the boundary. Wavelet filters sum to zero, so detail coefficients
% have zero mean and the estimator is their mean square. The scaling coefficients carry
% the mean of the series, which is not variance, so it is removed for that band.
wvar = zeros(1, level + 1);
for j = 1:level
	Lj = (2^j - 1) * (filtLen - 1) + 1;
	Mj = N - Lj + 1; % at least 1, by the cap on the level
	wvar(j) = mean(w(j, N - Mj + 1:end).^2);
end
Mj = N - ((2^level - 1) * (filtLen - 1) + 1) + 1;
wvar(level + 1) = var(w(level + 1, N - Mj + 1:end), 1);
numLevelsReturned = level; % number of detail levels

% ------------------------------------------------------------------------------
%% Return statistics
% ------------------------------------------------------------------------------
totalVar = sum(wvar);
out.scalingFrac = wvar(end) / totalVar; % variance fraction in the trend/scaling band

[~, out.domlevel] = max(wvar); % level (1..numLevelsReturned+1) carrying the most variance

% Scaling exponent: log2(variance) vs. level, across detail bands only
detailVar = wvar(1:numLevelsReturned);
validLevels = find(detailVar > 0);
if numel(validLevels) >= 2
	p = polyfit(validLevels, log2(detailVar(validLevels)), 1);
	out.decaySlope = p(1);
else
	out.decaySlope = NaN;
end

end
