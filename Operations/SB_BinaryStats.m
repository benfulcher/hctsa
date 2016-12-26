function out = SB_BinaryStats(y,binaryMethod)
% SB_BinaryStats    Statistics on a binary symbolization of the time series
%
% Binary symbolization of the time series is a symbolic string of 0s and 1s.
%
% Provides information about the coarse-grained behavior of the time series
%
%---INPUTS:
% y, the input time series
%
% binaryMethod, the symbolization rule:
%         (i) 'diff': by whether incremental differences of the time series are
%                      positive (1), or negative (0),
%          (ii) 'mean': by whether each point is above (1) or below the mean (0)
%          (iii) 'iqr': by whether the time series is within the interquartile range
%                      (1), or not (0).
%
%---OUTPUTS:
% Include the Shannon entropy of the string, the longest stretches of 0s
% or 1s, the mean length of consecutive 0s or 1s, and the spread of consecutive
% strings of 0s or 1s.

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
% Check inputs, set defaults:
%-------------------------------------------------------------------------------
if nargin < 2 || isempty(binaryMethod)
    binaryMethod = 'diff';
end

%-------------------------------------------------------------------------------
% Binarize the time series:
%-------------------------------------------------------------------------------
yBin = BF_Binarize(y,binaryMethod);

N = length(yBin); % length of signal - 1 (difference operation)

%-------------------------------------------------------------------------------
% Stationarity of binarized time series:
%-------------------------------------------------------------------------------
% (cf. SB_MotifTwo for basic stats on binarized time series)

% Stationarity:
out.pupstat2 = sum(yBin(floor(end/2)+1:end) == 1)/sum(yBin(1:floor(end/2)) == 1);

%-------------------------------------------------------------------------------
% Consecutive string of ones / zeros (normalized by length)
%-------------------------------------------------------------------------------
difffy = diff(find([1;yBin;1]));
stretch0 = difffy(difffy ~= 1) - 1;

difffy = diff(find([0;yBin;0] == 0));
stretch1 = difffy(difffy ~= 1) - 1;

%-------------------------------------------------------------------------------
% pstretches
%-------------------------------------------------------------------------------
% Number of different stretches as proportion of time series
out.pstretch1 = length(stretch1)/N;
% The following are trivially dependent on pstretch1:
% out.pstretch0 = length(stretch0)/N;
% out.pstretches = (length(stretch0)+length(stretch1))/N;

if isempty(stretch0) % all 1s (almost impossible to actually occur)
    out.longstretch0 = 0;
    out.meanstretch0 = 0;
    out.stdstretch0 = NaN;
else
    out.longstretch0 = max(stretch0); % longest consecutive stretch of zeros
    out.meanstretch0 = mean(stretch0); % mean stretch of zeros
    out.stdstretch0 = std(stretch0); % standard deviation of stretch lengths of consecutive zeros
end

if isempty(stretch1) % all 0s (almost impossible to actually occur)
    out.longstretch1 = 0;
    out.meanstretch1 = 0;
    out.stdstretch1 = NaN;
else
    out.longstretch1 = max(stretch1); % longest consecutive stretch of ones
    out.meanstretch1 = mean(stretch1);
    out.stdstretch1 = std(stretch1);
end

out.meanstretchdiff = out.meanstretch1 - out.meanstretch0;
out.stdstretchdiff = out.stdstretch1 - out.stdstretch0;

out.diff21stretch1 = mean(stretch1 == 2) - mean(stretch1 == 1);
out.diff21stretch0 = mean(stretch0 == 2) - mean(stretch0 == 1);

end
