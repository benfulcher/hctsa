function out = CO_autocorr(y,tau)
% INPUTS:
% y a scalar time series column vector
% If tau is a scalar, returns autocorrelation for y at that lag
% If tau is a vector, returns autocorrelations for y at that set of lags
% Added vector functionality for tau. Ben Fulcher 12/11/2009


N = length(y); % length of the time sries

if length(tau) == 1
    out = sum((y(1:N-tau) - mean(y(1:N-tau))).*(y(tau+1:N) ...
	            - mean(y(tau+1:N))))/N/std(y(1:N-tau))/std(y(tau+1:N));
else
    out = zeros(length(tau),1);
    for i = 1:length(tau)
        out(i) = sum((y(1:N-tau(i)) - mean(y(1:N-tau(i)))).*(y(tau(i)+1:N) ...
			        - mean(y(tau(i)+1:N))))/N/std(y(1:N-tau(i)))/std(y(tau(i)+1:N));
    end
end

end