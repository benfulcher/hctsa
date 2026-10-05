function out = DN_WithinP(x, p, meanOrMedian)
% DN_WithinP   Proportion of data points within a distance of the center of the distribution.
%
% Returns the proportion of data points that lie within p units of the center
% of the distribution. With 'mean', the center is the mean and the unit is the
% standard deviation. With 'median', the center is the median and the unit is
% the interquartile range divided by 1.35 (equal to the standard deviation for
% Gaussian data).
%
% ---INPUTS:
% x, the input data vector
% p, the number of units on each side of the center (default: 1)
% meanOrMedian, the center and unit to use (default: 'mean'):
%           'mean': the mean and standard deviation
%           'median': the median and iqr(x)/1.35
%
% ---OUTPUTS:
% a scalar: the proportion of data points within p units of the center.

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

% -------------------------------------------------------------------------------
% Check inputs:
% -------------------------------------------------------------------------------
if nargin < 2 || isempty(p)
	p = 1; % 1 std from mean
end
if nargin < 3 || isempty(meanOrMedian)
	meanOrMedian = 'mean';
end

% -------------------------------------------------------------------------------
% Compute the property:
% -------------------------------------------------------------------------------

N = length(x); % length of the time series

switch meanOrMedian
	case 'mean'
		mu = mean(x); % mean of the time series
		sig = std(x); % standard deviation of the time series

	case 'median'
		mu = median(x); % median of the time series
		sig = iqr(x) / 1.35; % rescaled interquartile range of the time series (equal
		% to standard deviation for Gaussian distribution)
	otherwise
		error('Unknown setting: ''%s''', meanOrMedian);
end

% The withinp statistic:
out = sum(x >= mu - p * sig & x <= mu + p * sig) / N;

end
