% CO_AutoCorr
% 
% Computes the autocorrelation of an input time series, y, at a time-lag, tau
% 
% INPUTS:
% y, a scalar time series column vector
% tau, the time-delay. If tau is a scalar, returns autocorrelation for y at that
%       lag. If tau is a vector, returns autocorrelations for y at that set of
%       lags.
%       
% Output is the autocorrelation at the given time-lag
% 

function out = CO_AutoCorr(y,tau)

% Check inputs:
if nargin < 2 || isempty(tau)
    tau = 1;
end
N = length(y); % length of the time sries

if length(tau) == 1 % output a single value at the given time-lag
    out = sum((y(1:N-tau) - mean(y(1:N-tau))).*(y(tau+1:N) ...
	            - mean(y(tau+1:N))))/N/std(y(1:N-tau))/std(y(tau+1:N));
else % output values over a range of time-lags
    out = zeros(length(tau),1);
    for i = 1:length(tau)
        out(i) = sum((y(1:N-tau(i)) - mean(y(1:N-tau(i)))).*(y(tau(i)+1:N) ...
			        - mean(y(tau(i)+1:N))))/N/std(y(1:N-tau(i)))/std(y(tau(i)+1:N));
    end
end

end