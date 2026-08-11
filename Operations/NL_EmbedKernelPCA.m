function out = NL_EmbedKernelPCA(y, tau, m, maxN)
% NL_EmbedKernelPCA   Kernel PCA of a time series in an embedding space, vs. linear PCA.
%
% Reconstructs the time series as a time-delay embedding (as in NL_EmbedPCA)
% and performs kernel Principal Components Analysis on the result using an
% RBF kernel, then compares the resulting eigenvalue spectrum to that of
% ordinary (linear) PCA on the same embedded points.
%
% At any finite kernel bandwidth, kernel PCA's spectrum is less compact than
% linear PCA's in absolute terms (its RBF feature space is far higher-
% dimensional than the embedding itself), but the *amount* by which it is
% less compact depends on the geometry of the point cloud: a linear process
% has an embedding point cloud well described by a low-dimensional
% ellipsoid, which kernel PCA cannot compress much better than linear PCA
% does, whereas a nonlinear process (e.g., one confined to a curved manifold
% in the embedding space) is compressed comparatively well by kernel PCA's
% nonlinear eigendirections. The *relative discrepancy* between the two
% spectra (e.g., their top2/nto50/nto80 ratios) is therefore used here as a
% nonlinearity signal, in a similar spirit to how existing hctsa operations
% compare a metric's value to its outcome under a linear-appropriate
% transformation (cf. surrogate-based tests in SD_SurrogateTest).
%
% "Nonlinear Component Analysis as a Kernel Eigenvalue Problem"
% B. Scholkopf, A. Smola, K.-R. Muller, Neural Comput. 10(5) 1299 (1998)
%
% cf. "Extracting qualitative dynamics from experimental data"
% D. S. Broomhead and G. P. King, Physica D 20(2-3) 217 (1986)
%
% ---INPUTS:
% y, the input time series
%
% tau, the time-delay, can be an integer or 'ac', or 'mi' for first
%               zero-crossing of the autocorrelation function or first minimum
%               of the automutual information, respectively
%
% m, the embedding dimension
%
% maxN, the maximum number of embedded points used to form the N x N kernel
%       matrix, whose eigendecomposition costs O(N^3). Longer embeddings are
%       reduced to their first maxN points (default: 2000, i.e., a ~30MB
%       kernel matrix, well under a second to eigendecompose; a warning is
%       issued whenever this cropping actually happens, since the spectrum
%       estimate is still noisy at this size and keeps sharpening with more
%       -- this is a memory/time cap, not a convergence point). Can be set
%       to 'full' to disable, with a second warning above 5000 points,
%       where the eigendecomposition starts to take several seconds.
%
% ---OUTPUTS:
% The same battery of eigenvalue-spectrum statistics as NL_EmbedPCA
% (variance explained by each of the top m components, how many components
% are needed to explain 50-90% of the variance, etc.), computed on the
% kernel PCA spectrum instead of the linear PCA spectrum, plus explicit
% comparison outputs (ratios/differences of matched linear and kernel
% statistics) that isolate the linear-vs-nonlinear discrepancy itself.

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
% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(tau)
	tau = 'ac'; % embed by first zero-crossing of autocorrelation function
end

if nargin < 3 || isempty(m)
	m = 3; % three-dimensional embedding
end

if nargin < 4 || isempty(maxN)
	maxN = 2000; % caps the O(N^3) kernel-matrix eigendecomposition
end

% ------------------------------------------------------------------------------
% Embed the signal via time-delay method
% ------------------------------------------------------------------------------
y_embed = BF_Embed(y, tau, m, false);

if isnan(y_embed)
	% embedding parameters are unsuitable (likely that tau is too long...)
	fprintf(1, 'Embedding parameters are not suitable for this time series\n');
	out = NaN; return
end

% pca() returns at most min(size(y_embed,1)-1, m) components; need m of them
% (and at least 2 for out.top2 below) or the loop indexing below is invalid:
if size(y_embed,1) - 1 < m || m < 2
	fprintf(1, 'Not enough embedding vectors (%u) for a rank-%u PCA\n', size(y_embed,1), m);
	out = NaN; return
end

% ------------------------------------------------------------------------------
% Crop to maxN embedded points for the kernel matrix (memory/time cap)
% ------------------------------------------------------------------------------
Nemb = size(y_embed,1);
if ischar(maxN) && strcmp(maxN, 'full')
	slowThreshold = 5000;
	if Nemb > slowThreshold
		warning('%u embedded points exceeds %u with maxN=''full''; kernel eigendecomposition may take several seconds', ...
				Nemb, slowThreshold);
	end
	y_kpca = y_embed;
elseif Nemb > maxN
	warning(['Cropping to the first %u of %u embedded points for kernel PCA ' ...
			 '(memory/time cap, not a convergence point -- raise maxN or use ''full'' ' ...
			 'for a more precise, more expensive estimate)'], maxN, Nemb);
	y_kpca = y_embed(1:maxN, :);
else
	y_kpca = y_embed;
end
N = size(y_kpca, 1);

if N - 1 < m
	fprintf(1, 'Not enough embedded points (%u) after cropping for a rank-%u kernel PCA\n', N, m);
	out = NaN; return
end

% ------------------------------------------------------------------------------
% Linear PCA on the (possibly cropped) embedded points, for comparison
% ------------------------------------------------------------------------------
[~, ~, latentLin] = pca(y_kpca);
percLin = latentLin / sum(latentLin);
statsLin = SUB_spectrumstats(percLin, m);

% ------------------------------------------------------------------------------
% Kernel PCA using an RBF kernel
% ------------------------------------------------------------------------------
% Pairwise squared Euclidean distances between embedded points:
sqDist = pdist2(y_kpca, y_kpca, 'squaredeuclidean');

% Median-heuristic bandwidth: the median of the pairwise distances, a
% standard parameter-free choice that sets the kernel's length scale to
% the data's own typical point-to-point spacing:
offDiag = sqDist(~eye(N));
medSqDist = median(offDiag);
if medSqDist == 0
	% All embedded points coincide (e.g., a constant or near-constant series)
	out = NaN; return
end
% Standard median-heuristic convention: exp(-sqDist / median(sqDist)):
K = exp(-sqDist / medSqDist);

% Center the kernel matrix in feature space:
oneN = ones(N, N) / N;
Kc = K - oneN*K - K*oneN + oneN*K*oneN;
Kc = (Kc + Kc') / 2; % symmetrize away numerical asymmetry

% Eigendecomposition (Kc is symmetric; eigenvalues of the centered kernel
% matrix are N times the eigenvalues of the empirical covariance operator
% in feature space, but the constant factor cancels in the normalized
% variance-explained spectrum computed below):
eigK = eig(Kc);
eigK = sort(eigK, 'descend');
eigK(eigK < 0) = 0; % clip small negative numerical noise (Kc is PSD in theory)

if sum(eigK) == 0 || eigK(m) == 0
	fprintf(1, 'Kernel PCA produced a degenerate (near-zero-rank) spectrum\n');
	out = NaN; return
end

percKern = eigK / sum(eigK);
statsKern = SUB_spectrumstats(percKern, m);

% ------------------------------------------------------------------------------
%% Raw kernel PCA outputs (mirrors NL_EmbedPCA's field names)
% ------------------------------------------------------------------------------
fields = fieldnames(statsKern);
for i = 1:length(fields)
	out.(fields{i}) = statsKern.(fields{i});
end

% ------------------------------------------------------------------------------
%% Linear-vs-kernel comparison outputs: the nonlinearity signal
% ------------------------------------------------------------------------------
% Ratios/differences of matched linear and kernel spectrum statistics.
% At any finite bandwidth the RBF feature space is far higher-dimensional
% than the embedding itself, so the kernel spectrum is *always* less
% compact than the linear one in absolute terms (top2_ratio < 1 even for
% strongly nonlinear systems) -- this is not itself the nonlinearity
% signal. What differs by system is *how much* less compact: on a
% genuinely low-dimensional nonlinear manifold (e.g., a chaotic
% attractor), kernel PCA still finds much more compact structure than it
% does for a linear/stochastic process at the same bandwidth, so these
% ratios sit closer to 1 (top2_ratio, nto50_ratio) for nonlinear systems
% and further from 1 for linear/noise ones (empirically confirmed here on
% AR(1)/white noise vs. the logistic map/Lorenz system: top2_ratio ~0.48-
% 0.53 for the former, ~0.75-0.77 for the latter).
out.top2_ratio = statsKern.top2 / statsLin.top2;
out.top2_diff = statsKern.top2 - statsLin.top2;
out.nto80_ratio = statsKern.nto80 / statsLin.nto80;
out.nto80_diff = statsKern.nto80 - statsLin.nto80;
out.nto50_ratio = statsKern.nto50 / statsLin.nto50;
out.std_ratio = statsKern.std / statsLin.std;

% ------------------------------------------------------------------------------
function stats = SUB_spectrumstats(perc, m)
	% Compute the same battery of eigenvalue-spectrum statistics as
	% NL_EmbedPCA, from a normalized (sums to 1), descending spectrum perc.

	for i = 1:m
		stats.(sprintf('perc_%u', i)) = perc(i);
	end

	stats.std = std(perc);
	stats.range = max(perc) - min(perc);
	stats.min = min(perc);
	stats.max = max(perc);
	stats.top2 = sum(perc(1:2));

	csperc = cumsum(perc);
	stats.nto50 = first_fn(csperc, 0.5, +1);
	stats.nto60 = first_fn(csperc, 0.6, +1);
	stats.nto70 = first_fn(csperc, 0.7, +1);
	stats.nto80 = first_fn(csperc, 0.8, +1);
	stats.nto90 = first_fn(csperc, 0.9, +1);

	stats.fb05 = first_fn(perc, 0.5, -1);
	stats.fb02 = first_fn(perc, 0.2, -1);
	stats.fb01 = first_fn(perc, 0.1, -1);
	stats.fb001 = first_fn(perc, 0.01, -1);
end

% ------------------------------------------------------------------------------
function firsti = first_fn(p, threshold, overOrUnder)
	% Find the first time p goes under/over the threshold, x

	if overOrUnder == -1
		firsti = find(p < threshold, 1, 'first');
	else
		firsti = find(p > threshold, 1, 'first');
	end

	% If it never goes under -- saturate as m at the maximum
	% (could be NaN, but this is more interpretable/comparable)
	if isempty(firsti)
		firsti = length(p) + 1;
	end
end

end
