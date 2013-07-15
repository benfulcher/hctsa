function out = ST_cumrange(y)
% Time series must be at least 50 points long
% Ben Fulcher, September 2009

doplot = 0; % plot outputs
N = length(y); % length of the time series
cums = zeros(N,1);

for i = 1:N
    cums(i) = range(y(1:i));
end
% cums=cums/range(y);

if doplot, figure('color','w'); plot(cums); end

fullr = range(y);

lunique = @(x) length(unique(x)); % use this alot -- how many unique entries

out.totnuq = lunique(cums);

% how many of the unique extrema are in first <proportion> of time series
cumtox = @(x) lunique(cums(1:floor(N*x)))/out.totnuq;
out.nuqp1 = cumtox(0.01);
out.nuqp10 = cumtox(0.1);
out.nuqp20 = cumtox(0.2);
out.nuqp50 = cumtox(0.5);
% out.nuqp1 = lunique(cums(1:floor(N*0.01)))/out.totnuq;
% out.nuqp10 = lunique(cums(1:floor(N*0.1)))/out.totnuq;
% out.nuqp20 = lunique(cums(1:floor(N*0.2)))/out.totnuq;
% out.nuqp50 = lunique(cums(1:floor(N*0.5)))/out.totnuq;

% how many unique extrema are in first <length> of time series
out.nuql10 = lunique(cums(1:10))/out.totnuq;
out.nuql50 = lunique(cums(1:50))/out.totnuq;

if N > 100
    out.nuql100 = lunique(cums(1:100))/out.totnuq;
else
    out.nuql100 = NaN;
end
    
if N > 1000
    out.nuql1000 = lunique(cums(1:1000))/out.totnuq;
else
    out.nuql1000 = NaN;
end


% (**2**) Actual proportion of full range captured at different points

out.p1 = cums(ceil(N*0.01))/fullr;
out.p10 = cums(ceil(N*0.1))/fullr;
out.p20 = cums(ceil(N*0.2))/fullr;
out.p50 = cums(ceil(N*0.5))/fullr;

out.l10 = cums(10)/fullr;
out.l50 = cums(50)/fullr;

if N > 100
    out.l100 = cums(100)/fullr;
else
    out.l100 = NaN;
end

if N > 1000,
    out.l1000 = cums(1000)/fullr;
else
    out.l1000 = NaN;
end

end