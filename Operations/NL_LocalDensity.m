function out = NL_LocalDensity(y, NNR, past, embedParams)
% NL_LocalDensity   How densely the delay-embedded trajectory is sampled around each of its points, and how that density changes along the orbit.
%
% Computes a k-nearest-neighbor estimate of the local probability density at each
% point of the time-delay embedding: density(i) = (k/Neff) / (V_m * r_k(i)^m),
% where r_k(i) is the distance from point i to its k-th (k = NNR) nearest
% neighbor (excluding temporally-close points within a Theiler window of "past"
% samples), m is the embedding dimension, V_m = pi^(m/2)/Gamma(m/2+1) is the volume
% of the unit m-ball, and Neff = N_embed - 2*past - 1 is the number of points that
% can be neighbors. The estimate is computed in units of the series' standard
% deviation, and its logarithm is analyzed:
%   log density(i) = log(k/Neff) - log(V_m) - m*log(r_k(i)/std(y)),
% which makes the outputs independent of the units of y (rescaling the series only
% shifts the log density, in a way that is removed by measuring distances in
% standard deviations), and of the number of points for a stationary process.
% Working with the log density rather than the density itself keeps the
% statistics well behaved (the density is heavy-tailed, with a few very dense
% points dominating its mean and standard deviation). To avoid infinite densities
% when there are repeated values (zero neighbor distance, as in quantized or
% held series), distances are smoothed as sqrt(r^2 + (0.01*median(r(r>0)))^2).
% This operation previously used TSTOOL's 'localdensity', which the original
% author noted was "very poorly documented in the TSTOOL package" -- its
% exact algorithm was never confirmed; the estimate is now computed natively
% (no toolbox dependency at all).
%
% The result is a series of log-density values in the time order of the embedded
% points, of length N - (m-1)*tau. The outputs describe its distribution and its
% serial dependence.
%
% ---INPUTS:
% y, the time series as a column vector
% NNR, number of nearest neighbours to compute (default: 3)
% past, Theiler window of time-correlated points to discard: {'ac', k} for k times
%       the first zero-crossing of the autocorrelation function, or a number of
%       samples (see BF_TheilerWindow; default: {'ac',1})
% embedParams, the embedding parameters, inputs to BF_Embed as {tau,m}, where
%              tau and m can be characters specifying a given automatic method
%              of determining tau ('ac', 'ac1e' or 'mi'; see BF_GetTau) and/or m
%              ('fnn') (see BF_Embed; default: {'ac','fnn'})
%
% ---OUTPUTS: statistics of the log local density series (output names retain 'den'):
% minden, maxden, iqrden, rangeden, stdden, meanden, medianden: minimum, maximum,
%       interquartile range, range, standard deviation, mean and median
% ac1den, ac2den, ac3den, ac4den, ac5den: autocorrelation at lags 1 to 5
% tauacden: the first zero-crossing of the autocorrelation function (with
%       interpolation)
% taumigaussden: the first minimum of the automutual information (the Gaussian
%       estimate, as CO_FirstMin(y,'mi-gaussian'), a monotonic function of the
%       autocorrelation)
%
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
%% Check inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(NNR)
	NNR = 3; % 3 nearest neighbours
end

if nargin < 3 || isempty(past)
	past = {'ac', 1};
end
past = BF_TheilerWindow(y, past);
if isnan(past) % the autocorrelation function never crosses zero
	warning('No autocorrelation zero-crossing to set the Theiler window')
	out = NaN; return
end

if nargin < 4 || isempty(embedParams)
	embedParams = {'ac', 'fnn'};
	fprintf(1, 'Using default embedding using autocorrelation and cao''s method.\n');
end

% ------------------------------------------------------------------------------
%% Embed the signal (native MATLAB matrix embedding, not a TSTOOL/TISEAN call)
% ------------------------------------------------------------------------------
Y = BF_Embed(y, embedParams{1}, embedParams{2}, false);

if isscalar(Y) && isnan(Y) % embedding failed
	warning('Embedding failed');
	out = NaN; return
end
[N_embed, m] = size(Y);

if N_embed <= NNR + 2 * past
	warning('Time series too short to do a local density estimate with these parameters');
	out = NaN; return
end

% ------------------------------------------------------------------------------
%% k-nearest-neighbor local density estimate
% ------------------------------------------------------------------------------
% Over-fetch candidate neighbors via a KD-tree, then discard any within the
% Theiler window (including the point itself); expand to a full (Theiler-
% window-excluding) search only for the rare point where that isn't enough:
kFetch = min(N_embed - 1, NNR + 2 * past + 5);
[idx, dist] = knnsearch(Y, Y, 'K', kFetch + 1);

dk = zeros(N_embed, 1); % distance to the NNR-th neighbor outside the Theiler window
for i = 1:N_embed
	validDists = dist(i, abs(idx(i, :) - i) > past);
	if length(validDists) < NNR
		allDists = sqrt(sum((Y - Y(i, :)).^2, 2));
		allDists(abs((1:N_embed)' - i) <= past) = Inf;
		validDists = sort(allDists);
	end
	dk(i) = validDists(NNR);
end

if ~any(dk > 0) % all neighbor distances are zero (e.g., a constant series)
	out = NaN; return
end

% Smooth the distances so that repeated values (zero distances) give a finite density:
d = sqrt(dk.^2 + (0.01 * median(dk(dk > 0)))^2);

% Log of the k-NN density estimate, with distances in units of the series' SD:
Neff = N_embed - 2 * past - 1; % number of points that can be neighbors of a given point
locden = log(NNR / Neff) - ((m / 2) * log(pi) - gammaln(m / 2 + 1)) - m * log(d / std(y));
% locden is a vector of length equal to the number of points in the
% embedding space (length of time series - (m-1)*tau), the log local density
% estimate at each point

out.minden = min(locden);
out.maxden = max(locden);
out.iqrden = iqr(locden);
out.rangeden = range(locden);
out.stdden = std(locden);
out.meanden = mean(locden);
out.medianden = median(locden);

F_acden = @(x) CO_AutoCorr(locden, x, 'Fourier'); % autocorrelation of locden for 1:5
for i = 1:5
	out.(sprintf('ac%uden', i)) = F_acden(i);
end

% Estimates of correlation length:
out.tauacden = CO_FirstCrossing(locden, 'ac', 0, 'continuous'); % first zero-crossing of autocorrelation function
out.taumigaussden = CO_FirstMin(locden, 'mi-gaussian'); % first minimum of the Gaussian automutual information function

end
