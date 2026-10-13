function yTransform = BF_RemovePoints(y, removeHow, p, removeOrSaturate, randomSeed)
% BF_RemovePoints   Remove or saturate a proportion of the points of a time series.
%
% Chooses a proportion, p, of the points of the (z-scored) series according to a rule,
% and either deletes them or clips their values. Removing deletes the chosen points and
% closes up the rest into a shorter series. Saturating keeps them in place but clips
% their values to the most extreme value among the points kept. Used by DN_RemovePoints
% (order-free statistics of the changed series) and CO_RemovePoints (autocorrelation
% statistics of the changed series).
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
% ---OUTPUTS:
% yTransform, the series after removing (a shorter series, in the original order) or
%       saturating (a series of the original length) the chosen points

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
N = length(y); % time-series length

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
	randomSeed = []; % default for BF_RandomSeed
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
		is = BF_Random(N, BF_RandomSeed(randomSeed), 'perm')'; % random ordering, reproducible
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

end
