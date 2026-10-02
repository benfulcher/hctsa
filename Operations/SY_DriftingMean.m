function out = SY_DriftingMean(y, segmentHow, l)
% SY_DriftingMean   Mean and variance in local time-series subsegments.
%
% Splits the time series into consecutive segments, computes the mean and variance
% in each segment, and compares the maximum and minimum segment means to the mean
% of the segment variances. A drifting mean makes the extreme segment means large
% relative to the within-segment variance. A final partial segment is dropped.
% Returns NaN if the segments are longer than the series.
%
% The idea is from a comp.soft-sys.matlab (MATLAB newsgroup) posting by Rune
% ("It seems to me that you are looking for a measure for a drifting mean. If so,
% this is what I would try: decide on a frame length N; split your signal in a
% number of frames of length N; compute the means of each frame; compute the
% variance for each frame; compare the ratio of maximum and minimum mean with the
% mean variance of the frames.")
%
% ---INPUTS:
% y, the input time series
%
% segmentHow, how to segment the series:
%       (i) 'fix': fixed-length segments (of length l)
%       (ii) 'num': a given number, l, of segments (default)
%
% l, either the length ('fix') or number ('num') of segments (default: 5 segments
%       for 'num', 200 samples for 'fix')
%
% ---OUTPUTS:
% max, the maximum segment mean divided by the mean of the segment variances
% min, the minimum segment mean divided by the mean of the segment variances
% mean, the mean of the segment means divided by the mean of the segment variances
% meanmaxmin, the average of max and min
% meanabsmaxmin, the average of the absolute values of max and min

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

N = length(y); % length of the input time series

% ------------------------------------------------------------------------------
%% Check inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(segmentHow)
	segmentHow = 'num'; % a specified number of segments
end
% (l's default must be set BEFORE the 'num' branch below converts it to a segment
% length -- that branch reads l, so with nargin < 3 it previously errored on an
% undefined variable rather than using the documented default.)
if nargin < 3 || isempty(l)
	switch segmentHow
		case 'num'
			l = 5; % 5 segments
		case 'fix'
			l = 200; % 200-sample segments
	end
end

if strcmp(segmentHow, 'num')
	l = floor(N / l);
elseif ~strcmp(segmentHow, 'fix')
	error('Unknown input setting ''%s''', segmentHow)
end

% ------------------------------------------------------------------------------
%% Check for short time series
% -------------------------------------------------------------------------------
if l == 0 || N < l % doesn't make sense to split into more windows than there are data points
	fprintf(1, 'Time Series (N = %u < l = %u) is too short for this operation\n', N, l);
	out = NaN;
	return
end

% -------------------------------------------------------------------------------
%% Get going
% -------------------------------------------------------------------------------
numFits = floor(N / l); % number of times l fits completely into N
z = zeros(l, numFits);
for i = 1:numFits
	z(:, i) = y((i - 1) * l + 1:i * l);
end
zm = mean(z);
zv = var(z);
meanVar = mean(zv);
maxMean = max(zm);
minMean = min(zm);
meanMean = mean(zm);

% -------------------------------------------------------------------------------
%% Output statistics
% -------------------------------------------------------------------------------

out.max = maxMean / meanVar;
out.min = minMean / meanVar;
out.mean = meanMean / meanVar;
out.meanmaxmin = (out.max + out.min) / 2;
out.meanabsmaxmin = (abs(out.max) + abs(out.min)) / 2;

end
