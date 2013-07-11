function glscf = CO_glscf(y,alpha,beta,tau)
% Calculates the generalized self-correlation function of a time series, as 
% introduced by Duarte Queiros and Moyano in Physica A, Vol. 383, pp. 10--15 (2007)
% The paper is titled "Yet on statistical properties of traded volume: 
% Correlation and mutual information at different value magnitudes"
% This function considers magnitude correlations
% Inputs alpha, beta are real and nonzero
% When alpha = beta estimates how values of the same order of magnitude are related in time
% When alpha ~= beta, estimates correlations between different magnitudes of the time series.
% Coded by Ben Fulcher, September 2009

%% Defaults
if nargin < 4 || isempty(tau)
    tau = 'tau';
end

% Set tau to first zero-crossing of the autocorrelation function with the input 'tau'
if strcmp(tau,'tau')
    tau = CO_fzcac(y);
end

y1 = abs(y(1:end-tau));
y2 = abs(y(1+tau:end));

glscf = (mean((y1.^alpha).*(y2.^beta)) - mean(y1.^alpha)*mean(y2.^beta)) / ...
 		(sqrt(mean(y1.^(2*alpha)) - mean(y1.^alpha)^2) * sqrt(mean(y2.^(2*beta)) - mean(y2.^beta)^2));


end