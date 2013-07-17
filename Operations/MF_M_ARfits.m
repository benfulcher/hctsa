function out = MF_M_ARfits(y,p)
% AR fit of order p
% Uses the arcov function from Matlab's Signal Processing Toolbox
% Ben Fulcher

N = length(y); % Length of input time series

if nargin < 2 || isempty(p)
    p = 2; % Fit AR(2) model by default
end

% Fit an AR model:
[a, e] = arcov(y,p);

out.e = e; % variance

% Output fitted parameters up to order, p (+1)
for i = 1:p+1
	eval(sprintf('out.a%u = a(%u);',i,i);
end

%% Residual analysis
y_est = filter([0, -a(2:end)],1,y);
err = y - y_est; % residuals

out.rms = sqrt(sum(err.^2)/N); % RMS error
out.mu = mean(err); % mean error
out.std = std(err); % std of error
out.AC1 = CO_autocorr(err,1); % autocorrelation of residuals at lag 1
out.AC2 = CO_autocorr(err,2); % autocorrelation of residuals at lag 2

end