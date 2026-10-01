function out = WL_Coeffs(y, wname, level)
% WL_Coeffs   How quickly the sorted magnitudes of a wavelet detail signal decay from their maximum.
%
% Performs a discrete wavelet decomposition (wavedec) of the time series using a
% given wavelet down to a given level, reconstructs the detail signal at that
% level alone (wrcoef, which has the same length N as the series), and sorts its
% magnitudes from largest to smallest. For each of a set of fractions p of the
% maximum, it returns the position in the sorted list, as a proportion of N, where
% the magnitudes first fall below p times the maximum. Small values indicate that
% a few large excursions dominate the detail at that scale.
%
% Uses Matlab's Wavelet Toolbox.
%
% ---INPUTS:
% y, the input time series
% wname, the wavelet name, e.g., 'db3' (see Wavelet Toolbox Documentation for
%        all options; default: 'db3')
% level, the level of wavelet decomposition (an integer, or 'max' for the maximum
%        level given by wmaxlev; default: 3). A level too large for the series
%        returns NaN.
%
% ---OUTPUTS:
% wb99m, wb90m, wb75m, wb50m, wb25m, wb10m, wb1m: the position (as a proportion of
%        the series length) in the sorted detail magnitudes at which they first
%        fall below 99%, 90%, 75%, 50%, 25%, 10%, 1% of their maximum ('where below
%        _ of maximum'). NaN if they never do.
%
% ---NOTES:
% The mean, maximum and median of the detail coefficients are not returned, since
% they duplicate WL_DWTCoeff's per-level stdd/noisestd fields (r > 0.98 on real
% data); the decay-shape profile above is this function's distinct contribution.

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
BF_CheckToolbox('wavelet_toolbox');

% ------------------------------------------------------------------------------
%% Check Inputs
N = length(y); % time-series length

if nargin < 2 || isempty(wname)
	wname = 'db3'; % default wavelet
end
if nargin < 3 || isempty(level)
	level = 3; % level of wavelet decomposition
end
if strcmp(level, 'max')
	level = wmaxlev(N, wname);
	if level == 0
		error('Cannot compute wavelet coefficients (short time series?)');
	end
end

if wmaxlev(N, wname) < level
	fprintf(1, 'Chosen level is too large for this wavelet on this signal\n');
	out = NaN;
	return
end

% ------------------------------------------------------------------------------
%% Perform a single-level wavelet decomposition
% (Recover a noisy signal by suppressing an approximation)
[c, l] = wavedec(y, level, wname);

% Reconstruct detail
det = wrcoef('d', c, l, wname, level); % detail this level

det_s = sort(abs(det), 'descend'); % sorted detail coefficient magnitudes

% plot(det_s);

% ------------------------------------------------------------------------------
%% Return statistics
% (mean/max/median of the detail coefficients are not returned here: they
% duplicate WL_DWTCoeff's per-level stdd/noisestd fields, which correlate
% r>0.98 with them on real data. WL_Coeffs' distinct contribution is the
% decay-shape profile below.)

% Decay rate stats ('where below _ maximum' = 'wb_m')
out.wb99m = findMyThreshold(0.99);
out.wb90m = findMyThreshold(0.90);
out.wb75m = findMyThreshold(0.75);
out.wb50m = findMyThreshold(0.50);
out.wb25m = findMyThreshold(0.25);
out.wb10m = findMyThreshold(0.10);
out.wb1m = findMyThreshold(0.01);

% ------------------------------------------------------------------------------
function propt = findMyThreshold(x)
	% where drops below proportion x of maximum
	propt = find(det_s < x * max(det_s), 1, 'first') / N;
	% (as a proportion of time-series length)
	if isempty(propt)
		propt = NaN;
	end
end
% ------------------------------------------------------------------------------

end
