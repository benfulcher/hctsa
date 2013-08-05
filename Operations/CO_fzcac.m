% CO_fzcac
% 
% Returns the first zero-crossing of the autocorrelation function.
% Uses CO_AutoCorr to calculate autocorrelations.
% 
% INPUTS:
% y, the input time series
% maxtau, a maximum time-delay to search up to
% 

function out = CO_fzcac(y,maxtau)
% Ben Fulcher, 2008

N = length(y); % the length of the time series

if nargin < 2 || isempty(maxtau)
    maxtau = N; % search up to a maximum of the length of the time series
    % maxtau = 400; % searches up to this maximum time lag
    % maxtau = min(maxtau,N); % searched up to the length of the time series if this is less than maxtau
end

% Calculate autocorrelation at increasing lags, until you find a negative one
for tau = 1:maxtau-1
    if CO_AutoCorr(y,tau) < 0
        out = tau; return
    end
end

% If haven't left yet, set output to maxtau
out = maxtau;

end