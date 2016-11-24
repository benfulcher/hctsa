function out = ST_LocalExtrema(y,lorf,n)
% ST_LocalExtrema   How local maximums and minimums vary across the time series.
%
% Finds maximums and minimums within given segments of the time series and
% analyses the results.
%
%---INPUTS:
% y, the input time series
%
% lorf, whether to use:
%     (i) 'l', windows of a given length (in which case the third input, n
%             specifies the length)
%     (ii) 'n', a specified number of windows to break the time series up into
%               (in which case the third input, n specifies this number)
%     (iii) 'tau', sets a window length equal to the correlation length of the
%                 time series, the first zero-crossing of the autocorrelation
%                 function.
%
% n, somehow specifies the window length given the setting of lorf above.

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

if nargin < 2 || isempty(lorf)
    lorf = 'l';
end
if nargin < 3 || isempty(n)
    switch lorf
       case 'l'
           n = 100; % 100-sample windows
       case 'n'
           n = 5; % 5 windows
    end
end

doPlot = 0; % plot outputs to a figure

N = length(y); % length of time series

% ------------------------------------------------------------------------------
%% Set window length
% ------------------------------------------------------------------------------
switch lorf
    case 'l'
        wl = n; % window length
    case 'n'
        wl = floor(N/n); % number of windows
    case 'tau'
        % this may not be a good idea!
        wl = CO_FirstZero(y,'ac');
    otherwise
        error('Unknown method ''%s''',lorf);
end

if (wl > N) || (wl <= 1)
    % ++BF 19/3/2010: this is not suitable if window length longer than ts
    fprintf(1,'The window length is longer than the time-series length!\n');
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Buffer the time series
% ------------------------------------------------------------------------------
y_buff = buffer(y,wl); % no overlap
% each *column* is a window of samples
if y_buff(end) == 0
    y_buff = y_buff(:,1:end-1); % remove last window if zero-padded
end
nw = size(y_buff,2); % number of windows

% ------------------------------------------------------------------------------
%% Find local extrema
% ------------------------------------------------------------------------------
locmax = max(y_buff); % summary of local maxima
locmin = min(y_buff); % summary of local minima
abslocmin = abs(locmin); % absoluate value of local minima
exti = find(abslocmin>locmax);
locext = locmax; locext(exti) = locmin(exti); % local extrema (furthest from mean; either maxs or mins)
abslocext = abs(locext); % the magnitude of the most extreme events in each window

if doPlot
    figure('color','w'); hold on;
    plot(locmax);
    plot(abslocext,'--g');
    plot(abslocmin,':r')
    plot(locext,'k');
end

% ------------------------------------------------------------------------------
%% Return outputs
% ------------------------------------------------------------------------------
out.meanrat = mean(locmax)/mean(abslocmin);
out.medianrat = median(locmax)/median(abslocmin);
out.minmax = min(locmax);
out.minabsmin = min(abslocmin);
out.minmaxonminabsmin = min(locmax)/min(abslocmin);
out.meanmax = mean(locmax);
out.meanabsmin = mean(abslocmin);
out.meanext = mean(locext);
out.medianmax = median(locmax);
out.medianabsmin = median(abslocmin);
out.medianext = median(locext);
out.stdmax = std(locmax);
out.stdmin = std(locmin);
out.stdext = std(locext);
out.zcext = ST_SimpleStats(locext,'zcross');
out.meanabsext = mean(abslocext);
out.medianabsext = median(abslocext);
out.diffmaxabsmin = sum(abs(locmax-abslocmin))/nw;
out.uord = sum(sign(locext))/nw; % whether extreme events are more up or down
out.maxmaxmed = max(locmax)/median(locmax);
out.minminmed = min(locmin)/median(locmin);
out.maxabsext = max(abslocext)/median(abslocext);


end
