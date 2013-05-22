function out = CO_fzcglscf(y,alpha,beta)
% Calculates the first minimum of the self-correlation function of a time series, as introduced by Duarte Queiros and Moyano in Physica A, Vol. 383, pp. 10--15 (2007)
% The paper is titled "Yet on statistical properties of traded volume: Correlation and mutual information at different value magnitudes"
% Keeps calculating until the function finds a minimum, and returns this lag (to a maximum of 400)
% Uses the function CO_glscf to calculate the generalized self-correlation function
% Ben Fulcher, 2009
% Ben Fulcher 2013--Now uses CO_glscf to calculate self-correlations

N = length(y); % the length of the time series

maxtau = 400; % searches up to this maximum time lag
maxtau = min(maxtau,N); % make sure no longer than the time series itself

glscfs = zeros(maxtau,1);

for i = 1:maxtau
	tau = i;
    % y1 = abs(y(1:end-tau));
    % y2 = abs(y(1+tau:end));
    
    glscfs(i) = CO_glscf(y,alpha,beta,tau);
    % glscfs(i) = (mean((y1.^alpha).*(y2.^beta)) - mean(y1.^alpha) * mean(y2.^beta)) / ...
    %                  (sqrt(mean(y1.^(2*alpha)) - mean(y1.^alpha)^2) ...
    %                          * sqrt(mean(y2.^(2*beta)) - mean(y2.^beta)^2));

	if i > 1 && glscfs(i)*glscfs(i-1) < 0
		% Draw a straight line between these two and look at where hits zero
		out = i - 1 + glscfs(i)/(glscfs(i)-glscfs(i-1));
		return
	end
end

out = maxtau; % if the function hasn't exited yet, set output to maxtau

end