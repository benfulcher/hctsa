function out = NL_localdensity(y, NNR, past, embedParams)
% NL_localdensity     Local density estimates in the time-delay embedding space
%
% Computes a standard k-nearest-neighbor local density estimate at each
% point of the time-delay embedding: density(i) is proportional to
% 1/r_NNR(i)^m, where r_NNR(i) is the distance from point i to its NNR-th
% nearest neighbor (excluding temporally-close points within a Theiler
% window of "past" samples) and m is the embedding dimension. This
% operation previously used TSTOOL's 'localdensity', which the original
% author noted was "very poorly documented in the TSTOOL package" -- its
% exact algorithm was never confirmed, only assumed to be some form of
% local density estimate in the embedding space, which is what's computed
% here natively (no toolbox dependency at all). The missing normalizing
% constant (relating 1/r^m to a true probability density) is the same for
% every point in a given call, so it cancels out of all of the relative
% statistics below (min/max/std/mean/median/autocorrelation).
%
% ---INPUTS:
%
% y, the time series as a column vector
%
% NNR, number of nearest neighbours to compute
%
% past, number of time-correlated points to discard (samples)
%
% embedParams, the embedding parameters, inputs to BF_Embed as {tau,m}, where
%               tau and m can be characters specifying a given automatic method
%               of determining tau and/or m (see BF_Embed).
%
% ---OUTPUTS: various statistics on the local density estimates at each point in
% the time-delay embedding, including the minimum and maximum values, the range,
% the standard deviation, mean, median, and autocorrelation.

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
	past = 40;
end

if nargin < 4 || isempty(embedParams)
	embedParams = {'ac', 'fnnmar'};
	fprintf(1, 'Using default embedding using autocorrelation and cao''s method.\n');
end

% ------------------------------------------------------------------------------
%% Embed the signal (native MATLAB matrix embedding, not a TSTOOL/TISEAN call)
% ------------------------------------------------------------------------------
Y = BF_Embed(y, embedParams{1}, embedParams{2}, 0);

if isscalar(Y) && isnan(Y) % embedding failed
	error('Embedding failed.')
end
[N_embed, m] = size(Y);

if N_embed <= NNR + 2 * past
	error('Time series too short to do a local density estimate with these parameters.')
end

% ------------------------------------------------------------------------------
%% k-nearest-neighbor local density estimate
% ------------------------------------------------------------------------------
% Over-fetch candidate neighbors via a KD-tree, then discard any within the
% Theiler window (including the point itself); expand to a full (Theiler-
% window-excluding) search only for the rare point where that isn't enough:
kFetch = min(N_embed - 1, NNR + 2 * past + 5);
[idx, dist] = knnsearch(Y, Y, 'K', kFetch + 1);

locden = zeros(N_embed, 1);
for i = 1:N_embed
	validDists = dist(i, abs(idx(i, :) - i) > past);
	if length(validDists) < NNR
		allDists = sqrt(sum((Y - Y(i, :)).^2, 2));
		allDists(abs((1:N_embed)' - i) <= past) = Inf;
		validDists = sort(allDists);
	end
	locden(i) = 1 / (validDists(NNR)^m);
end

if all(locden == 0) || any(~isfinite(locden))
	out = NaN; return
end
% locden is a vector of length equal to the number of points in the
% embedding space (length of time series - (m-1)*tau), the local density
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
out.taumiden = CO_FirstMin(locden, 'mi'); % first minimum of automutual information function

end
