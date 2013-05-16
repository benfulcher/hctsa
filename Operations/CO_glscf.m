function glscf = CO_glscf(y,alpha,beta,tau)
% Calculates the generalized self-correlation function, as introduced by
% Duarte Queiros and Moyano in Physica A, 2007
% Looks at magnitude correlations
% Inputs alpha, beta are real and nonzero
% When alpha=beta estimates how values of the same order of magnitude are related in time
% When alpha ~= beta, estimates correlations between different magnitudes of the time series.
% Coded by Ben Fulcher September 2009

if strcmp(tau,'tau'), tau = CO_fzcac(y); end

y1 = abs(y(1:end-tau));
y2 = abs(y(1+tau:end));

glscf = (mean((y1.^alpha).*(y2.^beta)) - mean(y1.^alpha)*mean(y2.^beta)) / ...
 		( sqrt(mean(y1.^(2*alpha)) - mean(y1.^alpha)^2) * sqrt(mean(y2.^(2*beta)) - mean(y2.^beta)^2) );


end