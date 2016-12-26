function out = ST_MomentCorr(x,windowLength,wOverlap,mom1,mom2,whatTransform)
% ST_MomentCorr   Correlations between simple statistics in local windows of a time series.
%
% Idea to implement by Nick S. Jones.
%
%---INPUTS:
% x, the input time series
%
% windowLength, the sliding window length (can be a fraction to specify a proportion of
%       the time-series length)
%
% wOverlap, the overlap between consecutive windows as a fraction of the window
%       length,
%
% mom1, mom2: the statistics to investigate correlations between (in each window):
%               (i) 'iqr': interquartile range
%               (ii) 'median': median
%               (iii) 'std': standard deviation (about the local mean)
%               (iv) 'mean': mean
%
% whatTransform: the pre-processing whatTransformormation to apply to the time series before
%         analyzing it:
%               (i) 'abs': takes absolute values of all data points
%               (ii) 'sqrt': takes the square root of absolute values of all
%                            data points
%               (iii) 'sq': takes the square of every data point
%               (iv) 'none': does no whatTransformormation

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

doPlot = 0; % plot outputs

N = length(x); % number of samples in the input signal

%-------------------------------------------------------------------------------
% Check inputs, set defaults
%-------------------------------------------------------------------------------

% Sliding window length (samples)
if nargin < 2 || isempty(windowLength)
    windowLength = 0.02; % 2% of the time-series length
end
if windowLength < 1
    windowLength = ceil(N*windowLength);
end

% sliding window overlap length
if nargin < 3 || isempty(wOverlap)
    wOverlap = 1/5;
end
if wOverlap < 1 % specify a fraction OF THE WINDOW LENGTH
    wOverlap = floor(windowLength*wOverlap);
end

if nargin < 4 || isempty(mom1)
    mom1 = 'mean';
end
if nargin < 5 || isempty(mom2)
    mom2 = 'std';
end

if nargin < 6 || isempty(whatTransform)
    whatTransform = 'none';
end

% ------------------------------------------------------------------------------
% Apply the specified whatTransformormation:
% ------------------------------------------------------------------------------
switch whatTransform
    case 'abs'
        x = abs(x);
    case 'sq'
        x = x.^2;
    case 'sqrt'
        x = sqrt(abs(x));
    case 'none'
        % x = x;
    otherwise
        error('Unknown tranformation ''%s''',whatTransform)
end

% ------------------------------------------------------------------------------
% Create the windows:
% ------------------------------------------------------------------------------
x_buff = buffer(x,windowLength,wOverlap);
numWindows = (N/(windowLength-wOverlap)); % number of windows

if size(x_buff,2) > numWindows
    % fprintf(1,'Should have %u columns but we have %u: removing last one',numWindows,size(x_buff,2))
    x_buff = x_buff(:,1:end-1);
end % lose last point

% ok, now we have the sliding window ('buffered') signal, x_buff
% first calculate the first moment in all the windows (each column is a
% 'window' of the signal

M1 = SUB_calcmemoments(x_buff,mom1);
M2 = SUB_calcmemoments(x_buff,mom2);

R = corrcoef(M1,M2);
out.R = R(2,1); % correlation coefficient
out.absR = abs(R(2,1)); % absolute value of correlation coefficient
out.density = range(M1)*range(M2)/N; % density of points in M1--M2 space
out.mi = IN_MutualInfo(M1,M2,'gaussian');
% out.mi = BF_MutualInformation(M1,M2,'range','range',floor(sqrt(N)));
% out.mi = BF_MutualInformation(M1,M2,[0,1],[0,1],floor(sqrt(N)));
% this is a poor choice of bin number -- M1 and M2 are not length N

if doPlot
    figure('color','w');
    plot(M1,M2,'.k');
end

% ------------------------------------------------------------------------------
function moms = SUB_calcmemoments(x_buff,momType)
    switch momType
        case 'mean'
            moms = mean(x_buff);
        case 'std'
            moms = std(x_buff);
        case 'median'
            moms = median(x_buff);
        case 'iqr'
            moms = iqr(x_buff);
        otherwise
            error('Unknown statistic ''%s''',momType)
    end
end
% ------------------------------------------------------------------------------

end
