function out = CO_NonlinearAutocorr(y, taus, doAbs)
% CO_NonlinearAutocorr   A custom nonlinear autocorrelation of a time series.
%
% Nonlinear autocorrelations are of the form <x_i x_{i-\tau_1} x_{i-\tau_2}...>
% (the mean of the product of the series with several delayed copies of itself).
% The usual two-point autocorrelations are <x_i.x_{i-\tau}>.
%
% Assumes that all the taus are much less than the length of the time
% series, N, so that the means can be approximated as the sample means and the
% standard deviations approximated as the sample standard deviations and so
% the z-scored time series can simply be used straight-up.
%
% ---INPUTS:
% y, the z-scored time series (an Nx1 vector)
% taus, a vector of the time delays as above:
%       [2] computes <x_i x_{i-2}>
%       [1,2] computes <x_i x_{i-1} x_{i-2}>
%       [1,1,3] computes <x_i x_{i-1}^2 x_{i-3}>
%       [0,0,1] computes <x_i^3 x_{i-1}>
% doAbs, [optional] a boolean: if true, takes the absolute value of the product
%        before taking the mean, which is useful for an odd number of factors
%        (default: true if length(taus) is even, otherwise false)
%
% ---OUTPUTS:
% out, a scalar: the mean of the product (or of its absolute value).
%
% ---NOTES:
% (*) For an odd number of factors (i.e., an even length of taus) the result will be
%     near zero for time-reversible processes, due to fluctuations about the mean,
%     even for highly correlated signals; hence the doAbs default.
%
% (*) doAbs = true is really a different operation that can't be compared with
%     the values obtained from taking doAbs = false (i.e., for odd lengths of taus).
%
% (*) It can be helpful to look at nonlinearAC at each iteration.

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
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
%% Check inputs & set defaults:
% ------------------------------------------------------------------------------
if nargin < 3 || isempty(doAbs) % use default settings for doAbs
	if rem(length(taus), 2) == 1
		% Odd number of time-lags
		doAbs = false;
	else
		% Even number of time-lags
		doAbs = true; % take abs, otherwise will be a very small number
	end
end
% -------------------------------------------------------------------------------

N = length(y); % time-series length
tMax = max(taus); % the maximum delay time

% Compute the autocorrelation sum iteratively
nonlinearAC = y(tMax + 1:N);
for i = 1:length(taus)
	nonlinearAC = nonlinearAC .* y(tMax - taus(i) + 1:N - taus(i));
end

% -------------------------------------------------------------------------------
% Compute output
if doAbs
	out = mean(abs(nonlinearAC));
else
	out = mean(nonlinearAC);
end

end
