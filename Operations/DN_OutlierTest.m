function out = DN_OutlierTest(y, p, justMe)
% DN_OutlierTest   How the mean and spread of the values change when the most extreme values are trimmed.
%
% Removes the p% highest and p% lowest values of the series (2p% in total; values
% below the p-th percentile and above the (100-p)-th percentile, strictly) and
% returns the mean of the remaining values, and their standard deviation divided
% by the standard deviation of the full series. Time order plays no part. For a
% z-scored series the full series has mean 0 and standard deviation 1, so these
% are the changes in the mean and in the spread when the tails are cut.
%
% ---INPUTS:
% y, the input data vector (should be z-scored)
% p, the percentage of values to remove beyond each of the upper and lower
%       percentiles (default: 2)
% justMe [opt], return a single number instead of a structure:
%       'mean': the mean of the middle portion of the data
%       'std': the std of the middle portion of the data, relative to that of the
%       full series
%
% ---OUTPUTS:
% mean, the mean of the middle (100 - 2p)% of the data
% std, the standard deviation of the middle (100 - 2p)% of the data, divided by the
%       standard deviation of the full series

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
% Check Inputs:
% -------------------------------------------------------------------------------

if nargin < 2 || isempty(p)
	p = 2; % by default, remove 2% of values from upper and lower percentiles
end
if nargin < 3
	justMe = ''; % return a structure with both the mean and std
end

% -------------------------------------------------------------------------------
% Get going:
% -------------------------------------------------------------------------------
% mean of the middle (100-2*p)% of the data
out.mean = mean(y(y > prctile(y, p) & y < prctile(y, 100 - p)));

% std of the middle (100-2*p)% ofthe data
out.std = std(y(y > prctile(y, p) & y < prctile(y, 100 - p))) / std(y); % [although std(y) should be 1]

% Output just a specified element of the output structure:
if ~isempty(justMe)
	switch justMe
		case 'mean'
			out = out.mean;
		case 'std'
			out = out.std;
		otherwise
			error('Unknown option ''%s''', justMe);
	end
end

end
