function out = ST_LocalExtrema(y,howToWindow,n)
% ST_LocalExtrema   How local maximums and minimums vary across the time series.
%
% Finds maximums and minimums within given segments of the time series and
% analyses the results.
%
%---INPUTS:
% y, the input time series
%
% howToWindow, whether to use:
%     (i) 'l', windows of a given length (in which case the third input, n
%             specifies the length)
%     (ii) 'n', a specified number of windows to break the time series up into
%               (in which case the third input, n specifies this number)
%     (iii) 'tau', sets a window length equal to the correlation length of the
%                 time series, the first zero-crossing of the autocorrelation
%                 function.
%
% n, somehow specifies the window length given the setting of howToWindow above.

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

% ------------------------------------------------------------------------------
% Check Inputs
% ------------------------------------------------------------------------------

if nargin < 2 || isempty(howToWindow)
    howToWindow = 'l';
end
if nargin < 3 || isempty(n)
    switch howToWindow
       case 'l'
           n = 100; % 100-sample windows
       case 'n'
           n = 5; % 5 windows
    end
end

doPlot = false; % plot outputs to a figure

N = length(y); % length of time series

% ------------------------------------------------------------------------------
%% Set window length
% ------------------------------------------------------------------------------
switch howToWindow
    case 'l'
        windowLength = n; % window length
    case 'n'
        windowLength = floor(N/n); % number of windows
    case 'tau'
        % this may not be a good idea!
        windowLength = CO_FirstCrossing(y,'ac',0,'discrete');
    otherwise
        error('Unknown method ''%s''',howToWindow);
end

if (windowLength > N) || (windowLength <= 1)
    % This feature is unsuitable if the window length exceeds ts
    fprintf(1,'The window length is longer than the time-series length!\n');
    out = NaN;
    return
end

% ------------------------------------------------------------------------------
%% Buffer the time series
% ------------------------------------------------------------------------------
y_buff = buffer(y,windowLength); % no overlap
% Each *column* is a window of samples:
if y_buff(end) == 0
    y_buff = y_buff(:,1:end-1); % remove last window if zero-padded
end
numWindows = size(y_buff,2); % number of windows

% ------------------------------------------------------------------------------
%% Find local extrema
% ------------------------------------------------------------------------------
locMax = max(y_buff); % summary of local maxima
locMin = min(y_buff); % summary of local minima
absLocMin = abs(locMin); % absolute value of local minima
exti = find(absLocMin > locMax);
locExt = locMax;
locExt(exti) = locMin(exti); % local extrema (furthest from mean; either maxs or mins)
absLocExt = abs(locExt); % the magnitude of the most extreme events in each window

if doPlot
    figure('color','w');
    hold('on');
    plot(locMax);
    plot(absLocExt,'--g');
    plot(absLocMin,':r')
    plot(locExt,'k');
end

% ------------------------------------------------------------------------------
%% Return outputs
% ------------------------------------------------------------------------------
out.meanrat = mean(locMax)/mean(absLocMin);
out.medianrat = median(locMax)/median(absLocMin);
out.minmax = min(locMax);
out.minabsmin = min(absLocMin);
out.minmaxonminabsmin = min(locMax)/min(absLocMin);
out.meanmax = mean(locMax);
out.meanabsmin = mean(absLocMin);
out.meanext = mean(locExt);
out.medianmax = median(locMax);
out.medianabsmin = median(absLocMin);
out.medianext = median(locExt);
out.stdmax = std(locMax);
out.stdmin = std(locMin);
out.stdext = std(locExt);
out.zcext = ST_SimpleStats(locExt,'zcross');
out.meanabsext = mean(absLocExt);
out.medianabsext = median(absLocExt);
out.diffmaxabsmin = sum(abs(locMax - absLocMin))/numWindows;
out.uord = sum(sign(locExt))/numWindows; % whether extreme events are more up or down
out.maxmaxmed = max(locMax)/median(locMax);
out.minminmed = min(locMin)/median(locMin);
out.maxabsext = max(absLocExt)/median(absLocExt);


end
