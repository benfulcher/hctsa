function out = CO_AutoCorrX2(y, taus, whatDirection)
% CO_AutoCorrX2   Asymmetric squared cross-correlation of a time series.
%
% Computes a lag-resolved generalization of the 'leverage effect' correlation
% used to detect asymmetric volatility feedback (e.g., Bouchaud, Matacz &
% Potters, Phys. Rev. Lett. 87, 228701 (2001)): instead of the usual
% autocorrelation <x_t.x_{t+tau}>, one term is squared:
%
%   forward:  <x_t.x_{t+tau}^2>   -- does the (signed) value now predict the
%                                     squared (unsigned/energy) value later?
%   backward: <x_t^2.x_{t+tau}>   -- does the squared value now predict the
%                                     (signed) value later?
%
% For time-reversible, linear processes these two are equal (up to sampling
% noise); a systematic difference between them is a signature of nonlinear,
% time-irreversible structure such as volatility clustering with a leverage
% (asymmetric) feedback.
%
% Note that CO_NonlinearAutocorr(x_z,[0,tau],false) and
% CO_NonlinearAutocorr(x_z,[tau,tau],false) already compute the backward and
% forward statistics (respectively) at a single lag; this function computes
% both directions efficiently across a whole vector of lags for building
% lag-profile ('shape') features, cf. CO_AutoCorrX2Shape.
%
% ---INPUTS:
% y, the input time series (should be z-scored: zero mean, unit variance)
% taus, a vector of (non-negative, integer) time lags to compute the
%       statistic at. tau = 0 reduces to skewness, mean(y.^3), for both
%       directions.
% whatDirection, 'forward' (default): <x_t.x_{t+tau}^2>
%                'backward':          <x_t^2.x_{t+tau}>
%
% ---OUTPUT:
% A vector (same length as taus) of the requested statistic at each lag.

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
if nargin < 2 || isempty(taus)
	taus = 1;
end
if nargin < 3 || isempty(whatDirection)
	whatDirection = 'forward';
end

y = y(:); % ensure column vector
N = length(y);

if max(taus) > N - 1
	warning('Time lag %u is too long for time-series length %u', max(taus), N)
end
if any(taus < 0)
	error('taus must be non-negative (this is an asymmetric statistic -- use whatDirection to choose the direction)');
end

% ------------------------------------------------------------------------------
%% Compute the statistic at each lag:
% ------------------------------------------------------------------------------
out = zeros(size(taus));
for i = 1:length(taus)
	tau = taus(i);
	yLater = y(1 + tau:N);   % x_{t+tau}
	yEarlier = y(1:N - tau); % x_t
	switch whatDirection
		case 'forward'
			% <x_t.x_{t+tau}^2>
			out(i) = mean(yEarlier .* yLater.^2);
		case 'backward'
			% <x_t^2.x_{t+tau}>
			out(i) = mean(yEarlier.^2 .* yLater);
		otherwise
			error('Unknown whatDirection ''%s'' (should be ''forward'' or ''backward'')', whatDirection);
	end
end

end
