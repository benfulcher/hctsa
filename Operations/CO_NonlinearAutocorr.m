function out = CO_NonlinearAutocorr(y,taus,doAbs)
% CO_NonlinearAutocorr      A custom nonlinear autocorrelation of a time series.
%
% Nonlinear autocorrelations are of the form:
% <x_i x_{i-\tau_1} x{i-\tau_2}...>
% The usual two-point autocorrelations are
% <x_i.x_{i-\tau}>
%
% Assumes that all the taus are much less than the length of the time
% series, N, so that the means can be approximated as the sample means and the
% standard deviations approximated as the sample standard deviations and so
% the z-scored time series can simply be used straight-up.
%
%---INPUTS:
% y  -- should be the z-scored time series (Nx1 vector)
% taus -- should be a vector of the time delays as above (mx1 vector)
%   e.g., [2] computes <x_i x_{i-2}>
%   e.g., [1,2] computes <x_i x_{i-1} x{i-2}>
%   e.g., [1,1,3] computes <x_i x_{i-1}^2 x{i-3}>
% doAbs [opt] -- a boolean (0,1) -- if one, takes an absolute value before
%                taking the final mean -- useful for an odd number of
%                contributions to the sum. Default is to do this for odd
%                numbers anyway, if not specified.
%
%---NOTES:
% (*) For odd numbers of regressions (i.e., even number length
%         taus vectors) the result will be near zero due to fluctuations
%         below the mean; even for highly-correlated signals. (doAbs)
%
% (*) doAbs = 1 is really a different operation that can't be compared with
%         the values obtained from taking doAbs = 0 (i.e., for odd lengths
%         of taus)
% (*) It can be helpful to look at nlac at each iteration.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
%% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 3 || isempty(doAbs) % use default settings for doAbs
    if rem(length(taus),2) == 1
        doAbs = 0;
    else
        % Even number of time-lags
        doAbs = 1; % take abs, otherwise will be a very small number
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
if doAbs
    out = mean(abs(nlac));
else
    out = mean(nlac);
end

end
