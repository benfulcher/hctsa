function out = CO_fzcac(y)
% Outputs the first zero-crossing of the autocorrelation function
% Uses CO_autocorr to calculate autocorrelations
% Very badly coded by Ben Fulcher, 2008
% Briefly tweaked by Ben Fulcher, May 2013

N = length(y); % the length of the time series

maxtau = 400; % searches up to this maximum time lag
maxtau = min(maxtau,N); % searched up to the length of the time series if this is less than maxtau

% Calculate autocorrelation at increasing lags, until you find a negative one
for tau = 1:maxtau-1
    if CO_autocorr(y,tau) < 0
        out = tau; return
    end
end
% If haven't left yet, set output to maxtau
out = maxtau;

end