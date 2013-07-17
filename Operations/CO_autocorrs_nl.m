% Computes autocorrelations of the input time series of the form
% <x_i x_{i-\tau_1} x{i-\tau_2}...>
% The usual two-point autocorrelations are
% <x_i.x_{i-\tau}>
% Assumes that all the taus are much less than N, the length of the time
% series so that the means can be approximated as the sample means and the
% standard deviations approximated as the sample standard deviations and so
% the zscored time series can simply be used straight up.
% Inputs:
% y  -- should be the z-scored time series (Nx1 vector)
% taus -- should be a vector of the time delays as above (mx1 vector)
%   e.g., [2] computes <x_i x_{i-2}>
%   e.g., [1,2] computes <x_i x_{i-1} x{i-2}>
%   e.g., [1,1,3] computes <x_i x_{i-1}^2 x{i-3}>
% mors [opt] -- take the mean or standard deviation of the final sum.
%               Standard deviation probably doesn't mean anything, but it
%               may well do! Default is of course to take the mean.
% doabs [opt] -- a boolean (0,1) -- if one, takes an absolute value before
%                taking the final mean -- useful for an odd number of
%                contributions to the sum. Default is to do this for odd
%                numbers anyway, if not specified.
%
% Note 1: for odd numbers of regressions (i.e., even number length
%         taus vectors) the result will be near zero due to fluctuations
%         below the mean; even for highly-correlated signals. (doabs)
% Note 2: doabs on really is a different metric that can't be compared with
%         the values obtained from taking doabs off (i.e., for odd lengths
%         of taus)  
% Note 3: It's nice to look at ely at each iteration.
% -- Ben Fulcher 8/6/09 --

function out = CO_autocorrs_nl(y,taus,mors,doabs)

%% Check Inputs:
if nargin < 3 || isempty(mors)
   mors = 'mean';
else
   if ~ismember(mors,{'mean','std'})
       error('Third input was ''%s'', should be either ''mean'' or ''std''.',mors);
   end
end
if nargin < 4 || isempty(doabs) % use default settings for doabs
    if rem(length(taus),2) == 1
        doabs = 1; % take abs, otherwise will be a very small number
    else
        doabs = 0; % not necessary to take absolute values
    end
end

N = length(y); % time-series length
tmax = max(taus); % the maximum delay time

% Compute the autocorrelation sum iteratively
ely = y(tmax+1:N);
for i = 1:length(taus)
    ely = ely.*y(tmax-taus(i)+1:N-taus(i));
end

% Compute output
switch mors
case 'mean'
    if doabs
        out = mean(ely);
    else
        out = mean(abs(ely));
    end
case 'std'
    if doabs
        out = std(ely);
    else
        out = std(abs(ely));
    end
end

end