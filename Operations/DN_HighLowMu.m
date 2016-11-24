function out = DN_HighLowMu(y)
% DN_HighLowMu      The highlowmu statistic.
%
% The highlowmu statistic is the ratio of the mean of the data that is above the
% (global) mean compared to the mean of the data that is below the global mean.
%
%---INPUTS:
% y, the input data vector

%---NOTES:
% Somehow measures the same information as SB_MotifTwo(y,'mean') -> u, i.e.,
% contains the same information as the proportion of the data that is above the
% mean. This indicates that you cannot independently control the proportion of
% data that is above the mean and the ratio of the means of the data above and
% below the mean. This is not immediately obvious...

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

mu = mean(y); % mean of data
mhi = mean(y(y > mu)); % mean of data above the mean
mlo = mean(y(y < mu)); % mean of data below the mean
out = (mhi-mu)/(mu-mlo); % ratio of the differences

end
