function out = CO_PartialAutoCorr(y, maxTau, whatMethod)
% CO_PartialAutoCorr   The partial autocorrelation of a time series.
%
% Computes the partial autocorrelation at lags 1 to maxTau: the correlation between
% y(t) and y(t-k) after removing the linear effect of the intermediate values (the
% last coefficient of an order-k autoregressive fit).
%
% The default ('burg') is the sequence of reflection coefficients of Burg's
% recursion: at each order, the coefficient minimizing the summed forward and backward
% prediction errors, which has a closed form and is bounded in [-1, 1]. It needs no
% matrix solve, so it is identical in any implementation, and it agrees closely with
% the ordinary-least-squares partial autocorrelation. The latter ('ols', MATLAB's
% parcorr) regresses y(t) on k lagged values; for smooth, nearly deterministic series
% its design matrix is numerically singular and the coefficient depends on the
% linear solver.
%
% ---INPUTS:
% y, a scalar time series column vector
% maxTau, the maximum time delay; returns lags up to this maximum (default 10)
% whatMethod, the method used to compute it: 'burg' (the default; Burg recursion) or
%               'ols' (ordinary least squares, parcorr)
%
% ---OUTPUTS:
% pac_1, pac_2, ..., pac_<maxTau>, the partial autocorrelation at lags 1, 2, ...,
%       maxTau (pac_1 to pac_20 for maxTau = 20).
%
% ---NOTES:
% For 'burg', each order k uses the N-k forward errors f(t) and the matching
% backward errors b(t-1) of the order k-1 fit: the reflection coefficient is
% 2*sum(f*b)/(sum(f^2) + sum(b^2)), the partial autocorrelation at lag k is its
% value, and both error series are then updated with it. A constant series gives
% zeros.
% The 'ols' method requires the Econometrics Toolbox (parcorr).

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
%% Check inputs and set defaults:
% ------------------------------------------------------------------------------
if nargin < 2
	% Use a maximum lag of 10 by default
	maxTau = 10;
end

if nargin < 3 || isempty(whatMethod)
	% Burg recursion by default
	whatMethod = 'burg';
end

% -------------------------------------------------------------------------------
%% Initial checks on maxTau
% -------------------------------------------------------------------------------
N = length(y); % time-series length

if maxTau <= 0
	error('Negative time lags not applicable')
end

% ------------------------------------------------------------------------------
%% Do the computation
% ------------------------------------------------------------------------------
switch whatMethod
case 'burg'
	nLags = min(maxTau, N - 1);
	f = y(:) - mean(y); % forward prediction errors (order 0: the series itself)
	b = f; % backward prediction errors
	% Prediction errors below 1e-8 of the series' energy are rounding-level noise; they
	% are not allowed to produce partial autocorrelations of up to +-1 (exactly
	% predictable series, such as a sinusoid, give zeros beyond their order)
	tiny = max(2e-8 * (f' * f), realmin);
	pacf = zeros(maxTau + 1, 1);
	pacf(1) = 1;
	for k = 1:nLags
		ff = f(k + 1:N);
		bb = b(k:N - 1);
		refl = 2 * (ff' * bb) / (ff' * ff + bb' * bb + tiny);
		f(k + 1:N) = ff - refl * bb;
		b(k + 1:N) = bb - refl * ff;
		pacf(k + 1) = refl;
	end
	pacf(nLags + 2:end) = NaN; % lags beyond the series length are undefined
case 'ols'
	pacf = parcorr(y, 'NumLags', maxTau, 'Method', 'ols');
otherwise
	error('Unknown method ''%s'' (use ''burg'' or ''ols'')', whatMethod);
end

% Zero lag is the first entry in the PACF (and should always be 1)
for i = 1:maxTau
	out.(sprintf('pac_%u', i)) = pacf(i + 1);
end

end
