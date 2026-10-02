function out = DN_RemovePoints(y, removeHow, p, removeOrSaturate, randomSeed)
% DN_RemovePoints   How the distribution of a time series changes when a set of points is removed or clipped.
%
% A proportion, p, of the points of the (z-scored) series are removed, or
% saturated, according to a rule (see BF_RemovePoints), and order-free statistics
% of the changed series are computed: its mean, median and standard deviation, the
% change in its skewness from that of the original series, and the ratio of its
% kurtosis to that of the original series. Removing
% deletes the chosen points and closes up the rest into a shorter series. Saturating
% keeps them in place but clips their values to the most extreme value among the
% points kept. The autocorrelation statistics of the same transformation are in
% CO_RemovePoints.
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
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed
%       (only relevant for removeHow = 'random', which is otherwise
%       irreproducible run to run; no registered feature uses it)
%
% ---OUTPUTS: statistics of the changed series, relative to the original:
% mean, median, std, the mean, median and standard deviation of the changed
%       series (not ratios; the z-scored original has mean 0 and std 1)
% skewnessdiff, the skewness of the changed series minus the skewness of the original
% kurtosisrat, the ratio of the kurtosis of the changed series to that of the
%       original
%
% ---NOTES:
% A similar idea is implemented in DN_OutlierInclude.
% The change in skewness is a difference rather than a ratio because the skewness of
% the original series can be near 0 (symmetric distributions), where a ratio is
% unstable (this output was previously the ratio, skewnessrat). The kurtosis
% is at least 1 for any series, so its ratio is always well defined.
% The autocorrelation outputs of this function (fzcacrat, ac1diff, ac2diff, ac3diff,
% sumabsacfdiff) moved to CO_RemovePoints.

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

% -------------------------------------------------------------------------------
%% Compute output statistics
% -------------------------------------------------------------------------------
out.mean = mean(yTransform);
out.median = median(yTransform);
out.std = std(yTransform);

% Requires Statistics Toolbox:
out.skewnessdiff = skewness(yTransform) - skewness(y);
out.kurtosisrat = kurtosis(yTransform) / kurtosis(y);

end
