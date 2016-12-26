function out = SY_StatAv(y,whatType,n)
% SY_StatAv     Simple mean-stationarity metric, StatAv.
%
% The StatAv measure divides the time series into non-overlapping subsegments,
% calculates the mean in each of these segments and returns the standard deviation
% of this set of means.
%
% cf. "Heart rate control in normal and aborted-SIDS infants", S. M. Pincus et al.
% Am J. Physiol. Regul. Integr. Comp. Physiol. 264(3) R638 (1993)
%
%---INPUTS:
%
% y, the input time series
%
% whatType, the type of StatAv to perform:
%           (i) 'seg': divide the time series into n segments
%           (ii) 'len': divide the time series into segments of length n
%
% n, either the number of subsegments ('seg') or their length ('len')

% Might be nicer to use the 'buffer' function for this...?
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

% ------------------------------------------------------------------------------
% Check Inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(whatType)
    whatType = 'seg'; % divide into n segments by default
end

if nargin < 3 || isempty(n)
    n = 5; % use 5 segments
end

N = length(y); % Time-series length

% ------------------------------------------------------------------------------
% Compute the means in local time-series segments
% ------------------------------------------------------------------------------

switch whatType
case 'seg'
    % divide time series into n segments
    M = zeros(n,1);
    p = floor(N/n);% lose the last N mod n data points

    for j = 1:n
        M(j) = mean(y(p*(j-1)+1:p*j));
    end

case 'len'
    if N > 2*n
        pn = floor(N/n);
        M = zeros(pn,1);
        for j = 1:pn
            M(j) = mean(y((j-1)*n+1:j*n));
        end
    else
        fprintf(1,'This time series (N = %u) is too short for StatAv(%s,''%u'')\n',N,whatType,n);
        out = NaN; return
    end

otherwise
    error('Error evaluating StatAv of type ''%s'', please select either ''seg'' or ''len''',whatType)

end

% ------------------------------------------------------------------------------
% Compute the statistic
% ------------------------------------------------------------------------------

s = std(y); % should be 1 (for a z-scored time-series input)
sdav = std(M);
out = sdav/s;

end
