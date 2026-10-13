function out = CO_RemovePoints(y, removeHow, p, removeOrSaturate, randomSeed)
% CO_RemovePoints   How the autocorrelation of a time series changes when a set of points is removed or clipped.
%
% A proportion, p, of the points of the (z-scored) series are removed, or
% saturated, according to a rule (see BF_RemovePoints), and the autocorrelation
% structure is compared before and after the change. Removing deletes the chosen
% points and closes up the rest into a shorter series, which splices together points
% that were not neighbors. Saturating keeps them in place but clips their values to
% the most extreme value among the points kept. The order-free statistics of the same
% transformation (mean, median, standard deviation, skewness and kurtosis) are in
% DN_RemovePoints.
%
% ---INPUTS:
% y, the input time series (should be z-scored)
% removeHow, how to choose the points to remove:
%       'absclose': those closest to the mean
%       'absfar': those furthest from the mean
%       'min': the lowest values
%       'max': the highest values
%       'random': at random
%       Default: 'absfar'.
% p, the proportion of points to remove (default: 0.1)
% removeOrSaturate, whether to remove the points ('remove', the default) or to
%       saturate their values ('saturate'; not possible with 'absclose' or
%       'random')
% randomSeed, the seed of the random ordering (see BF_RandomSeed; default: 0), which
%       comes from the portable generator BF_Random (only relevant for
%       removeHow = 'random'; no registered feature uses it)
%
% ---OUTPUTS: statistics of the changed series, relative to the original:
% fzcacrat, the ratio of the first zero-crossing of the autocorrelation function
%       (changed to original)
% ac1diff, ac2diff, ac3diff, the absolute differences in the autocorrelation at
%       lags 1, 2 and 3
% sumabsacfdiff, the sum over lags 1 to 8 of the absolute differences in the
%       autocorrelation
%
% ---NOTES:
% Only the first zero-crossing is a ratio: its original value is an interpolated lag
% of at least 0.5, so the ratio is always well defined. The autocorrelation outputs
% are differences, because the original autocorrelation can be near 0, where a ratio
% is unstable (the ratios ac1rat, ac2rat and ac3rat, redundant with the differences,
% have been removed).
% This function and DN_RemovePoints were split from a single function that returned
% both sets of outputs.

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
%% Check inputs
% ------------------------------------------------------------------------------
if nargin < 2
	removeHow = []; % defaults set in BF_RemovePoints
end
if nargin < 3
	p = [];
end
if nargin < 4
	removeOrSaturate = [];
end
if nargin < 5
	randomSeed = [];
end

% ------------------------------------------------------------------------------
%% Remove or saturate the chosen points
% ------------------------------------------------------------------------------
yTransform = BF_RemovePoints(y, removeHow, p, removeOrSaturate, randomSeed);

% Compute some autocorrelation properties:
acf_y = SUB_acf(y, 8);
acf_yTransform = SUB_acf(yTransform, 8);

% -------------------------------------------------------------------------------
%% Compute output statistics
% -------------------------------------------------------------------------------

% Two main comparison functions:
f_absDiff = @(x1, x2) abs(x1 - x2); % ignores the sign
f_ratio = @(x1, x2) x1 / x2;

out.fzcacrat = f_ratio(CO_FirstCrossing(yTransform, 'ac', 0, 'continuous'), ...
					   CO_FirstCrossing(y, 'ac', 0, 'continuous'));

out.ac1diff = f_absDiff(acf_yTransform(1), acf_y(1));

out.ac2diff = f_absDiff(acf_yTransform(2), acf_y(2));

out.ac3diff = f_absDiff(acf_yTransform(3), acf_y(3));

out.sumabsacfdiff = sum(abs(acf_yTransform - acf_y));

% -------------------------------------------------------------------------------
function acf = SUB_acf(x, n)
	% computes autocorrelation of the input sequence, x, up to a maximum time
	% lag, n
	acf = CO_AutoCorr(x, 1:n, 'Fourier');
end

end
