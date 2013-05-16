function out = DT_isseas(y)
% tests to see if seasonal by a simple sin model fit test
% Input: y: time series to test
% Output: 1 (seasonal), 0 (non-seasonal)
% Ben Fulcher 9/7/09

%% Foreplay
if size(y,2)>size(y,1)
    y=y';
end
y = zscore(y);% so that all time series are normalized to the same amplitude
N = length(y); % length of time series
r = [1:N]'; % range on which to fit

%% Fit
[cfun,gof] = fit(r,y,'sin1'); % fits form: a1*sin(b1*x+c1)

% Condition 1: fit is ok
th_fit=0.3; % r2>th_fit
% Condition 2: amplitude is not too small
th_ampl=0.5; % a1>th_ampl
% gof.rsquare
% cfun.a1
if gof.rsquare>th_fit && abs(cfun.a1>th_ampl)
    out=1;
else
    out=0;
end

end