function out = DN_pleft(y,th)
% DN_pleft      Distance from the mean at which a given proportion of data are more distant.
%
% Measures the maximum distance from the mean at which a given fixed proportion,
% p, of the time-series data points are further.
% Normalizes by the standard deviation of the time series
% (could generalize to separate positive and negative deviations in future)
% Uses the quantile function from Matlab's Statistics Toolbox
%
%---INPUTS:
% y, the input data vector
% th, the proportion of data further than p from the mean
%           (output p, normalized by standard deviation)

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

if nargin < 2 || isempty(th)
    th = 0.1; % default
end

p = quantile(abs(y-mean(y)),1-th);

% A proportion, th, of the data lie further than p from the mean
out = p/std(y);

end
