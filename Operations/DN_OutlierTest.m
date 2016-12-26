function out = DN_OutlierTest(y,p,justMe)
% DN_OutlierTest    How distributional statistics depend on distributional outliers.
%
% Removes the p% of highest and lowest values in the time series (i.e., 2*p%
% removed from the time series in total) and returns the ratio of either the
% mean or the standard deviation of the time series, before and after this
% transformation.
%
%---INPUTS:
% y, the input data vector
% p, the percentage of values to remove beyond upper and lower percentiles
% justMe [opt], just returns a number:
%               (i) 'mean' -- returns the mean of the middle portion of the data
%               (ii) 'std' -- returns the std of the middle portion of the data

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
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

%-------------------------------------------------------------------------------
% Check Inputs:
%-------------------------------------------------------------------------------

if nargin < 2 || isempty(p)
    p = 2; % by default, remove 2% of values from upper and lower percentiles
end
if nargin < 3
    justMe = ''; % return a structure with both the mean and std
end

%-------------------------------------------------------------------------------
% Get going:
%-------------------------------------------------------------------------------
% mean of the middle (100-2*p)% of the data
out.mean = mean(y(y > prctile(y,p) & y < prctile(y,100-p)));

% std of the middle (100-2*p)% ofthe data
out.std = std(y(y > prctile(y,p) & y < prctile(y,100-p))) / std(y); % [although std(y) should be 1]

% Output just a specified element of the output structure:
if ~isempty(justMe)
    switch justMe
    case 'mean'
        out = out.mean;
    case 'std'
        out = out.std;
    otherwise
        error('Unknown option ''%s''',justMe);
    end
end

end
