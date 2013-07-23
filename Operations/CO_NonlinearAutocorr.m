% CO_NonlinearAutocorr
% 
% Computes autocorrelations of the input time series of the form
% <x_i x_{i-\tau_1} x{i-\tau_2}...>
% The usual two-point autocorrelations are
% <x_i.x_{i-\tau}>
% 
% Assumes that all the taus are much less than the length of the time
% series, N, so that the means can be approximated as the sample means and the
% standard deviations approximated as the sample standard deviations and so
% the z-scored time series can simply be used straight-up.
% 
% % INPUTS:
% y  -- should be the z-scored time series (Nx1 vector)
% taus -- should be a vector of the time delays as above (mx1 vector)
%   e.g., [2] computes <x_i x_{i-2}>
%   e.g., [1,2] computes <x_i x_{i-1} x{i-2}>
%   e.g., [1,1,3] computes <x_i x_{i-1}^2 x{i-3}>
% doabs [opt] -- a boolean (0,1) -- if one, takes an absolute value before
%                taking the final mean -- useful for an odd number of
%                contributions to the sum. Default is to do this for odd
%                numbers anyway, if not specified.
%
% Note: for odd numbers of regressions (i.e., even number length
%         taus vectors) the result will be near zero due to fluctuations
%         below the mean; even for highly-correlated signals. (doabs)
% Note: doabs is really a different metric that can't be compared with
%         the values obtained from taking doabs off (i.e., for odd lengths
%         of taus)
% Note: It can be helpful to look at nlac at each iteration.

function out = CO_NonlinearAutocorr(y,taus,doabs)
% Ben Fulcher, 8/6/09

%% Check Inputs:
if nargin < 3 || isempty(doabs) % use default settings for doabs
    if rem(length(taus),2) == 1
        doabs = 0; 
    else
        % Even number of time-lags
        doabs = 1; % take abs, otherwise will be a very small number
    end
end

N = length(y); % time-series length
tmax = max(taus); % the maximum delay time

% Compute the autocorrelation sum iteratively
nlac = y(tmax+1:N);
for i = 1:length(taus)
    nlac = nlac.*y(tmax-taus(i)+1:N-taus(i));
end

% Compute output
if doabs
    out = mean(abs(nlac));
else
    out = mean(nlac);
end


end