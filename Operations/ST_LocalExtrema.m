% ST_LocalExtrema
% 
% Finds maximums and minimums within given segments of the time series and
% analyses the results. Outputs quantify how local maximums and minimums vary
% across the time series.
% 
% INPUTS:
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
% 

function out = ST_LocalExtrema(y,lorf,n)
% Ben Fulcher, August 2008

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

doplot = 0; % plot outputs to a figure

N = length(y); % length of time series

%% set window length
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

%% Buffer the time series
y_buff = buffer(y,wl); % no overlap
% each *column* is a window of samples
if y_buff(end) == 0
    y_buff = y_buff(:,1:end-1); % remove last window if zero-padded
end
nw = size(y_buff,2); % number of windows

%% Find local extrema
locmax = max(y_buff); % summary of local maxima
locmin = min(y_buff); % summary of local minima
abslocmin = abs(locmin); % absoluate value of local minima
exti = find(abslocmin>locmax);
locext = locmax; locext(exti) = locmin(exti); % local extrema (furthest from mean; either maxs or mins)
abslocext = abs(locext); % the magnitude of the most extreme events in each window

if doplot
    figure('color','w'); hold on;
    plot(locmax);
    plot(abslocext,'--g');
    plot(abslocmin,':r')
    plot(locext,'k');
end

%% Return outputs
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