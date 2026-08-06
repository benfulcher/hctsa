function out = WL_modwtvar(y, wname, level)
% WL_modwtvar   Multiscale variance decomposition via the maximal overlap DWT.
%
% Decomposes the time series into octave-scale bands using the maximal
% overlap discrete wavelet transform (MODWT) and computes summary statistics
% on how variance is distributed across scales. Unlike the standard
% (decimated) DWT used elsewhere in this codebase, the MODWT is shift-invariant
% and its associated variance estimator (modwtvar) is unbiased and accounts
% for boundary-affected coefficients at each level -- the standard approach
% for a scale-wise variance decomposition (Percival & Walden).
%
% ---INPUTS:
% y, the input time series
%
% wname, the mother wavelet, e.g., 'db3', 'sym2' (see Wavelet Toolbox
%           Documentation)
%
% level, the level of wavelet decomposition (can be set to 'max' for the
%               maximum level supported by the series length,
%               floor(log2(N)))
%
% ---OUTPUTS:
% scalingFrac, the fraction of variance in the lowest-frequency (scaling/trend)
%               band -- the part of the signal not resolved by any detail level
% domlevel, the level (1 = highest frequency, ..., level+1 = scaling/trend
%               band) that carries the most variance
% decaySlope, the slope of log2(variance) vs. level across the detail bands --
%               a wavelet-based scaling exponent, analogous to a Hurst
%               estimate but using the MODWT's unbiased, boundary-corrected
%               variance rather than an ad hoc regression on raw coefficients

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

maxLevelAllowed = floor(log2(N));
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
w = modwt(y, wname, level);
wvar = modwtvar(w, wname); % nominally length level+1: levels 1..level (detail), then scaling/trend
% modwtvar silently drops levels whose unbiased estimate has no
% boundary-unaffected coefficients left (short series relative to the
% wavelet filter length), so wvar can be shorter than level+1 -- always
% work from its actual length rather than the requested level.
numLevelsReturned = numel(wvar) - 1; % number of detail levels actually returned

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
