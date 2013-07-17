function out = CO_acfremove(y,howtorem,p)
% Analyze how the autocorrelation function changes as points are removed from a time series
% Remove a proportion p points from the full time series using a specified removal algorithm, howtorem
% Input time series, y, should be z-scored
% Ben Fulcher, September 2009

%% Preliminaries
N = length(y); % time-series length
doplot = 0; % plot output

%% Check inputs
if nargin < 2 || isempty(howtorem)
    howtorem = 'absfar'; % default
end
if nargin < 3 || isempty(p)
    p = 0.1; % 10%
end

switch howtorem
    case 'absclose' % remove a proportion p of points closest to the mean
        [~, is] = sort(abs(y),'descend');
    case 'absfar' % remove a proportion p of points furthest from the mean
        [~, is] = sort(abs(y),'ascend');
    case 'min'
        [~, is] = sort(y,'ascend');
    case 'max'
        [~, is] = sort(y,'descend');
    case 'random'
        is = randperm(N);
    otherwise
        error('Unknwon method ''%s''',howtorem);
end

rkeep = sort(is(1:round(N*(1-p))),'ascend');
y_trim = y(rkeep);

if doplot
    figure('color','w')
    hold off
    plot(y,'ok');
    hold on;
    plot(rkeep,y_trim,'.r')
    hold off;
    hist(y_trim,50)
end

acf_y = SUB_acf(y,8);
acf_y_trim = SUB_acf(y_trim,8);

if doplot
    figure('color','w')
    hold off;
    plot(acf_y,':b'); hold on;
    plot(acf_y_trim,':r');
end


%% Compute outputs
out.fzcacrat = CO_fzcac(y_trim)/CO_fzcac(y);
out.ac2rat = acf_y_trim(2)/acf_y(2); % includes the sign
out.ac2diff = abs(acf_y_trim(2)-acf_y(2));
out.ac3rat = acf_y_trim(3)/acf_y(3); % includes the sign
out.ac3diff = abs(acf_y_trim(3)-acf_y(3));
out.sumabsacfdiff = sum(abs(acf_y_trim-acf_y));
out.mean = mean(y_trim);
out.median = median(y_trim);
out.std = std(y_trim);
out.skewnessrat = skewness(y_trim)/skewness(y); % Statistics Toolbox
out.kurtosisrat = kurtosis(y_trim)/kurtosis(y); % Statistics Toolbox


function acf = SUB_acf(x,n)
    acf = zeros(n,1);
    for i = 1:n
        acf(i) = CO_autocorr(x,i-1);
    end
end

end