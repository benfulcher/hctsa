% DT_isseas
% 
% A simple test of seasonality by fitting a 'sin1' model to the time series
% using fit function from the Curve Fitting Toolbox. The output is binary: 1 if
% the goodness of fit, R^2, exceeds 0.3 and the amplitude of the fitted periodic
% component exceeds 0.5, and 0 otherwise.
% 
% INPUTS:
% y, the input time series
% 
% OUTPUT: 1 (seasonal), 0 (non-seasonal)

function out = DT_isseas(y)
% Ben Fulcher, 9/7/2009

%% Preliminaries
% Make sure the input time series, y, is a column vector
if size(y,2) > size(y,1)
    y = y';
end
N = length(y); % length of input time series
r = (1:N)'; % range over which to fit

%% Fit a sinusoidal model
[cfun, gof] = fit(r,y,'sin1'); % fits the following form: a1*sin(b1*x+c1)

%% Two conditions for determining whether time series contains periodicities:
% Condition 1: fit is ok
th_fit = 0.3; % r2>th_fit

% Condition 2: amplitude is not too small
th_ampl = 0.5; % a1>th_ampl

if gof.rsquare > th_fit && abs(cfun.a1 > th_ampl)
    out = 1; % test thinks the time series has strong periodicities
else
    out = 0; % test thinks the time series doesn't have any strong periodicities
end

end