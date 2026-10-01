function out = DN_RemovePoints(y, removeHow, p, removeOrSaturate, randomSeed)
% DN_RemovePoints   How time-series properties change when a set of points is removed or clipped.
%
% A proportion, p, of the points of the (z-scored) series are removed, or
% saturated, according to a rule, and statistics are computed before and after
% the change. Removing deletes the chosen points and closes up the rest into a
% shorter series. Saturating keeps them in place but clips their values to the
% most extreme value among the points kept.
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
% fzcacrat, the ratio of the first zero-crossing of the autocorrelation function
%       (changed to original)
% ac1rat, ac2rat, ac3rat, the ratios of the autocorrelation at lags 1, 2 and 3
%       (changed to original; the sign is kept)
% ac1diff, ac2diff, ac3diff, the absolute differences in the autocorrelation at
%       lags 1, 2 and 3
% sumabsacfdiff, the sum over lags 1 to 8 of the absolute differences in the
%       autocorrelation
% mean, median, std, the mean, median and standard deviation of the changed
%       series (not ratios; the z-scored original has mean 0 and std 1)
% skewnessrat, kurtosisrat, the ratios of the skewness and of the kurtosis
%       (changed to original)
%
% ---NOTES:
% A similar idea is implemented in DN_OutlierInclude.
% skewnessrat divides by the skewness of the original series, which is near 0 for
% symmetric distributions, so it is unstable for such series.

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
%% Preliminaries
% ------------------------------------------------------------------------------
N = length(y); % time-series length
doPlot = false; % plot output

% ------------------------------------------------------------------------------
%% Check inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(removeHow)
	removeHow = 'absfar'; % default
end
if nargin < 3 || isempty(p)
	p = 0.1; % 10%
end
if nargin < 4 || isempty(removeOrSaturate)
	removeOrSaturate = 'remove';
end

if nargin < 5
	randomSeed = []; % default for BF_ResetSeed
end

if ~BF_iszscored(y)
	warning('The input time series should be z-scored')
end

% ------------------------------------------------------------------------------
% Sort time-series values on different criteria, ordered by those to be *kept*
switch removeHow
	case 'absclose'
		% Remove a proportion p of points closest to the mean
		[~, is] = sort(abs(y), 'descend');
	case 'absfar'
		% Remove/saturate a proportion p of points furthest from the mean
		[~, is] = sort(abs(y), 'ascend');
	case 'min'
		% Remove/saturate a proportion p of points with the lowest values
		[~, is] = sort(y, 'descend');
	case 'max'
		% Remove/saturate a proportion p of points with the highest values
		[~, is] = sort(y, 'ascend');
	case 'random'
		BF_ResetSeed(randomSeed); % (for reproducibility of the random ordering)
		is = randperm(N);
	otherwise
		error('Unknown method ''%s''', removeHow);
end

% Indices of points to *keep*:
rKeep = sort(is(1:round(N * (1 - p))), 'ascend');

% Indices of points to *transform*:
rTransform = setxor(1:N, rKeep);

% -------------------------------------------------------------------------------
% Do the removing/saturating to convert y -> yTransform
switch removeOrSaturate
	case 'remove'
		% Remove the targeted points:
		yTransform = y(rKeep);

	case 'saturate'
		% Saturate out the targeted points:
		switch removeHow
			case 'max'
				yTransform = y;
				yTransform(rTransform) = max(y(rKeep));
			case 'min'
				yTransform = y;
				yTransform(rTransform) = min(y(rKeep));
			case 'absfar'
				yTransform = y;
				yTransform(yTransform > max(y(rKeep))) = max(y(rKeep));
				yTransform(yTransform < min(y(rKeep))) = min(y(rKeep));
			otherwise
				error('Cannot ''saturate'' when using ''%s'' method', removeHow)
		end
	otherwise
		error('Unknown removeOrSaturate option: ''%s''', removeOrSaturate);
end

% -------------------------------------------------------------------------------
% SIMPLE PLOT:
if doPlot
	figure('color', 'w')
	hold('off')
	plot(y, 'ok');
	hold('on');
	plot(rKeep, yTransform, '.r')
	hold('off');
	histogram(yTransform, 50)
end

% Compute some autocorrelation properties:
acf_y = SUB_acf(y, 8);
acf_yTransform = SUB_acf(yTransform, 8);

if doPlot
	figure('color', 'w')
	hold('off');
	plot(acf_y, ':b');
	hold('on');
	plot(acf_yTransform, ':r');
end

% -------------------------------------------------------------------------------
%% Compute output statistics
% -------------------------------------------------------------------------------

% Two main comparison functions:
f_absDiff = @(x1, x2) abs(x1 - x2); % ignores the sign
f_ratio = @(x1, x2) x1 / x2; % includes the sign

out.fzcacrat = f_ratio(CO_FirstCrossing(yTransform, 'ac', 0, 'continuous'), ...
					   CO_FirstCrossing(y, 'ac', 0, 'continuous'));

out.ac1rat = f_ratio(acf_yTransform(1), acf_y(1));
out.ac1diff = f_absDiff(acf_yTransform(1), acf_y(1));

out.ac2rat = f_ratio(acf_yTransform(2), acf_y(2));
out.ac2diff = f_absDiff(acf_yTransform(2), acf_y(2));

out.ac3rat = f_ratio(acf_yTransform(3), acf_y(3));
out.ac3diff = f_absDiff(acf_yTransform(3), acf_y(3));

out.sumabsacfdiff = sum(abs(acf_yTransform - acf_y));
out.mean = mean(yTransform);
out.median = median(yTransform);
out.std = std(yTransform);

% Requires Statistics Toolbox:
out.skewnessrat = skewness(yTransform) / skewness(y);
out.kurtosisrat = kurtosis(yTransform) / kurtosis(y);

% -------------------------------------------------------------------------------
function acf = SUB_acf(x, n)
	% computes autocorrelation of the input sequence, x, up to a maximum time
	% lag, n
	acf = CO_AutoCorr(x, 1:n, 'Fourier');
	% acf = zeros(n,1);
	% for i = 1:n
	%     acf(i) = CO_AutoCorr(x,i,'Fourier');
	% end
end

end
