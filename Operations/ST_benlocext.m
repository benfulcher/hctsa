function out = ST_benlocext(y,lorf,n)
% Finds local extrema in a given window of the time series and then does
% stuff with them :)
% Ben Fulcher August 2008

% lorf='l';n=100;
% y: the time series
% lorf: either 'l': length of window, 'n': number of windows
% , or 'tau': the autocorrelation length.
% n: the relevant number from lorf

N = length(y); % length of time series

%% set window length
switch lorf
    case 'l'
        wl = n; % window length
    case 'n'
        wl = floor(N/n); % number of windows
    case 'tau'
        % this may not be a good idea!
        wl = CO_fzcac(y);
end

if wl > N || wl <= 1
    % ++BF 19/3/2010: this is not suitable if window length longer than ts
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

% hold on; plot(locmax);  plot(abslocext,'--g'); plot(abslocmin,':r')
% plot(locext,'k');

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
out.zcext = ST_bs(locext,'zcross');
out.meanabsext = mean(abslocext);
out.medianabsext = median(abslocext);
out.diffmaxabsmin = sum(abs(locmax-abslocmin))/nw;
out.uord = sum(sign(locext))/nw; % whether extreme events are more up or down
out.maxmaxmed = max(locmax)/median(locmax);
out.minminmed = min(locmin)/median(locmin);
out.maxabsext = max(abslocext)/median(abslocext);


end