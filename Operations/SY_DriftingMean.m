function out = SY_DriftingMean(y,howl,l)
% SY_DriftingMean   Mean and variance in local time-series subsegments.
%
% Splits the time series into segments, computes the mean and variance in each
% segment and compares the maximum and minimum mean to the mean variance.
%
% This function implements an idea found in the Matlab Central forum:
% http://www.mathworks.de/matlabcentral/newsreader/view_thread/136539
%
% >> It seems to me that you are looking for a measure for a drifting mean.
% >> If so, this is what I would try:
% >>
% >> - Decide on a frame length N
% >> - Split your signal in a number of frames of length N
% >> - Compute the means of each frame
% >> - Compute the variance for each frame
% >> - Compare the ratio of maximum and minimum mean
% >>   with the mean variance of the frames.
% >>
% >> Rune
%
%---INPUTS:
% y, the input time series
%
% howl, (i) 'fix': fixed-length segments (of length l)
%       (ii) 'num': a given number, l, of segments
%
% l, either the length ('fix') or number of segments ('num')

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

N = length(y); % length of the input time series

% ------------------------------------------------------------------------------
% Check inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(howl)
    howl = 'num'; % a specified number of segments
end
if strcmp(howl,'num')
    l = floor(N/l);
elseif ~strcmp(howl,'fix')
    error('Unknown input setting ''%s''',howl)
end

if nargin < 3 || isempty(l)
    switch howl
    case 'num'
        l = 5; % 5 segments
    case 'fix'
        l = 200; % 200-sample segments
    end
end

% ------------------------------------------------------------------------------

% ++BF 19/3/2010
if N < l % doesn't make sense to split into more windows than there are data points
    fprintf(1,'Time Series (N = %u < l = %u) is too short for this operation\n',N,l);
    out = NaN; return
end

% Get going
nfits = floor(N/l);
z = zeros(l,nfits);
for i = 1:nfits; % number of times l fits completely into N
    z(:,i) = y((i-1)*l+1:i*l);
end
zm = mean(z);
zv = var(z);
meanvar = mean(zv);
maxmean = max(zm);
minmean = min(zm);
meanmean = mean(zm);

out.max = maxmean/meanvar;
out.min = minmean/meanvar;
out.mean = meanmean/meanvar;
out.meanmaxmin = (out.max+out.min)/2;
out.meanabsmaxmin = (abs(out.max)+abs(out.min))/2;

end
