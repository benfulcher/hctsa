function out = ST_cumrange(y)
% Time series must be at least 50 points long
% Ben Fulcher approx. September 2009.

N = length(y);
cums = zeros(N,1);

for i = 1:N
    cums(i) = range(y(1:i));
end
% cums=cums/range(y);

% plot(cums)

fullr = range(y);

out.totnuq = length(unique(cums));

% how many of the unique extrema are in first <proportion> of time series
out.nuqp1 = length(unique(cums(1:floor(N*0.01))))/out.totnuq;
out.nuqp10 = length(unique(cums(1:floor(N*0.1))))/out.totnuq;
out.nuqp20 = length(unique(cums(1:floor(N*0.2))))/out.totnuq;
out.nuqp50 = length(unique(cums(1:floor(N*0.5))))/out.totnuq;

% how many unique extrema are in first <length> of time series
out.nuql10 = length(unique(cums(1:10)))/out.totnuq;
out.nuql50 = length(unique(cums(1:50)))/out.totnuq;

if N > 100
    out.nuql100 = length(unique(cums(1:100)))/out.totnuq;
else
    out.nuql100 = NaN;
end
    
if N > 1000
    out.nuql1000 = length(unique(cums(1:1000)))/out.totnuq;
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