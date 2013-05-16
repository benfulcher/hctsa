function out = MF_M_ARfits(y,p)
% AR fit of order p
% Ben Fulcher

[a e] = arcov(y,p);

out.e = e; % variance

for i = 1:p+1
	eval(['out.a' num2str(i) '=a(i);']);
end

%% Residual analysis
y_est = filter([0 -a(2:end)],1,y);
N = length(y);
err = y-y_est;

out.rms = sqrt(sum(err.^2)/N); % rms error
out.mu = mean(err); % mean error
out.std = std(err); % std of error
out.AC1 = CO_autocorr(err,1); % autocorrelation of residuals at lag 1
out.AC2 = CO_autocorr(err,2); % autocorrelation of residuals at lag 2


end