function out = NL_EmbedCluster(y, tau, m, maxN)
% NL_EmbedCluster   Whether the time-delay embedding of the series forms separate clusters of points.
%
% Reconstructs the time series as a time-delay embedding (as in NL_EmbedPCA,
% NL_EmbedKernelPCA) and asks whether the resulting point cloud is made of
% distinct groups. A dynamical process whose trajectory visits distinct regions
% of phase space (e.g., alternating between two attractor states, or a system
% with intermittent bursts) leaves a multi-modal point cloud in the embedding; a
% process with a single smooth (e.g., unimodal-stochastic) attractor does not.
% This is a distinct signal from marginal-distribution multi-modality (e.g., DN_
% kurtosis/bimodality stats on the raw values), since two states can overlap
% entirely in amplitude yet still separate cleanly once lagged coordinates are
% added, and from regime-switching detected by MF_hmm_Fit / MF_hmm_CompareNStates,
% which cluster points in raw-amplitude (not lagged/embedded) space.
%
% All statistics are computed on the whitened cloud (zero mean, identity
% covariance), so they do not depend on its scale or orientation. sep_pp,
% sep_mahal and valley measure the separation of a two-group split of the cloud;
% they are always defined, even when two groups are not favored over one, and grow
% with the separation of genuine clusters. nongauss measures how far the cloud is
% from a single Gaussian, whatever the reason.
%
% ---INPUTS:
% y, the input time series
% tau, the time-delay, can be an integer, 'ac' for the first zero-crossing of the
%      autocorrelation function, 'ac1e' for the floor of its first 1/e crossing, or
%      'mi' for the smaller of the first minimum of the Kraskov automutual
%      information and the 'ac1e' delay (see BF_GetTau; default: 'ac')
% m, the embedding dimension (default: 2)
% maxN, the maximum number of embedded points used. Defaults to 'full' (no
%       cropping). Set to a number to reduce longer embeddings to their first maxN
%       points instead (a warning is issued whenever that cropping happens, since
%       the estimates keep sharpening with more points: it is a memory/time cap,
%       not a convergence point).
%
% ---OUTPUTS:
% sep_pp, projection-pursuit two-means separation: the largest fraction of the
%       variance of the whitened point cloud that lies between the two groups of
%       an optimal two-group split along a direction, over a fixed grid of
%       directions (see NOTES). About 2/pi (0.64) for Gaussian data, and 1 for
%       two distinct point masses.
% sep_mahal, log1p of the Mahalanobis separation between the two component
%       means of a 2-component Gaussian mixture fit, using their pooled
%       covariance (log-compressed to tame the heavy tail from near-singular
%       covariance on near-deterministic embeddings). The fit starts from
%       deterministic hard splits of the points (see NOTES).
% valley, the depth of the density valley at the sep_pp split: one minus the
%       kernel density of the projected points at the cut, relative to the
%       smaller of the densities at the two groups' median points, floored at 0.
%       0 when the projected density has no dip at the cut (e.g., any unimodal
%       cloud), and 1 when the two groups are separated by an empty gap.
% nongauss, the negentropy of the whitened cloud, in nats: the entropy of a
%       Gaussian with the same covariance minus the cloud's entropy, estimated
%       from nearest-neighbor distances (see NOTES). Near 0 for Gaussian data, and
%       larger the more structured (clustered, curved, skewed or heavy-tailed)
%       the cloud.
%
% ---NOTES:
% The cloud is whitened along its principal axes (each eigenvector signed so
% that its largest-magnitude element is positive).
%
% sep_pp: the whitened cloud is projected onto a fixed grid of directions: 90
% equally spaced angles in [0, pi) for m = 2; 300 Fibonacci-lattice points on the
% upper hemisphere for m = 3; the single axis for m = 1; and 300 points of a
% generalized golden-ratio (low-discrepancy) sequence mapped to the sphere for
% m > 3. For each direction the exact optimal one-dimensional two-means split is
% found (scanning every cut point between distinct sorted values), and its
% between-group sum of squares as a fraction of the total is recorded; sep_pp is
% the maximum over directions.
%
% sep_mahal: the 2-component mixture is fitted by EM from the hard split at the
% optimal one-dimensional two-means threshold along the first principal
% component, and again from the split found for sep_pp. Each fit starts from
% the two groups' proportions, means and regularized covariances (a
% regularization of 1e-6 times the mean variance of the embedded coordinates),
% and iterates until the mean log-likelihood per point changes by less than 1e-6
% (at most 500 iterations); the fit of higher likelihood is kept. If neither
% split can start a fit (e.g., a group with no more than m points), splits at the
% median of the same two projections are used. sep_mahal is NaN only if no
% 2-component fit is possible.
%
% valley: the points are projected onto the sep_pp direction (unit variance) and
% their density is estimated with a Gaussian kernel of bandwidth 0.9*N^(-1/5),
% Silverman's rule of thumb for unit variance. The density is evaluated at the
% cut and at each group's lower median point (a data point, so that both
% densities are positive).
%
% nongauss: the entropy of the whitened cloud is estimated with the
% Kozachenko-Leonenko nearest-neighbor estimator (k = 5 neighbors), with each
% nearest-neighbor distance floored at 0.01 (in whitened units), so that repeated
% values (point masses) give a finite, saturating value rather than a divergent
% one. The Gaussian reference entropy is that of the identity covariance,
% (m/2)*log(2*pi*e).
%
% A curved-but-unimodal manifold (e.g., the ring traced out by a periodic signal
% in a 2-d embedding) has a bimodal projection along any direction, and so also
% raises sep_pp, sep_mahal and valley, despite having no distinct dynamical
% states; it also raises nongauss, which measures departure from Gaussianity in
% general (e.g., a skewed marginal distribution also raises it), not separation
% specifically.

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
	m = 2; % two-dimensional embedding by default
end

if nargin < 4 || isempty(maxN)
	maxN = 'full'; % no cropping by default
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

% Need enough points, relative to m, for well-posed two-group fits and density
% estimates:
Nembed = size(y_embed, 1);
if Nembed < 80 * m
	fprintf(1, 'Not enough embedded points (%u) for stable %u-dimensional cluster statistics\n', ...
			Nembed, m);
	out = NaN; return
end

% A constant (or near-constant) embedded point cloud cannot be usefully clustered:
if all(range(y_embed) < 1e-10)
	out = NaN; return
end

% ------------------------------------------------------------------------------
% Crop to maxN embedded points (memory/time cap)
% ------------------------------------------------------------------------------
if ischar(maxN) && strcmp(maxN, 'full')
	y_cloud = y_embed;
else
	if Nembed > maxN
		warning(['Cropping to the first %u of %u embedded points ' ...
				 '(memory/time cap, not a convergence point -- raise maxN or use ''full'' ' ...
				 'for a more precise, more expensive estimate)'], maxN, Nembed);
		y_cloud = y_embed(1:maxN, :);
	else
		y_cloud = y_embed;
	end
end

% ------------------------------------------------------------------------------
% Whitened cloud (zero mean, identity covariance) along its principal axes
% ------------------------------------------------------------------------------
centered = y_cloud - mean(y_cloud, 1);
[evecs, evals] = eig(cov(centered));
[evals, ix] = sort(diag(evals), 'ascend');
evecs = evecs(:, ix);
% (fix the arbitrary sign of each eigenvector: its largest-magnitude element is positive)
[~, iMax] = max(abs(evecs), [], 1);
evecs = evecs .* sign(evecs(sub2ind([m, m], iMax, 1:m)));
evals = max(evals, 1e-12 * max(evals)); % (a rank-deficient cloud)
white = (centered * evecs) ./ sqrt(evals');
pc1 = white(:, end); % (the largest-variance axis; a rescaling does not change a split)

% ------------------------------------------------------------------------------
% sep_pp: best optimal two-means split over a fixed grid of directions
% ------------------------------------------------------------------------------
directions = SplitDirections(m);
numDirections = size(directions, 1);
r2 = zeros(numDirections, 1);
thr = zeros(numDirections, 1);
for lo = 1:50:numDirections % (chunked to bound memory for long series)
	hi = min(lo + 49, numDirections);
	[r2(lo:hi), thr(lo:hi)] = OptimalTwoMeansSplits(white * directions(lo:hi, :)');
end
[out.sep_pp, iBest] = max(r2);
ppProj = white * directions(iBest, :)';

% ------------------------------------------------------------------------------
% sep_mahal: 2-component mixture from deterministic hard-split starts (along the
% sep_pp direction and along PC1)
% ------------------------------------------------------------------------------
% Regularize covariance estimates proportionally to the data's own scale
% (y is z-scored, but the embedded coordinates can still have any variance):
regVal = 1e-6 * mean(var(y_cloud));
[~, pc1Thr] = OptimalTwoMeansSplits(pc1);
gm2 = [];
for attempt = 1:2
	if attempt == 1
		starts = {ppProj > thr(iBest), pc1 > pc1Thr};
	else % (fallback: splits at the median)
		starts = {ppProj > median(ppProj), pc1 > median(pc1)};
	end
	bestLL = -Inf;
	for j = 1:2
		gmj = GMM2FromLabels(y_cloud, starts{j}, regVal);
		if ~isempty(gmj) && gmj.meanLogLik > bestLL
			gm2 = gmj;
			bestLL = gmj.meanLogLik;
		end
	end
	if ~isempty(gm2)
		break
	end
end
if isempty(gm2)
	out.sep_mahal = NaN;
else
	pooledCov = (gm2.Sigma(:, :, 1) + gm2.Sigma(:, :, 2)) / 2;
	dMu = gm2.mu(1, :) - gm2.mu(2, :);
	% log1p-compressed: the raw Mahalanobis separation is numerically unstable for
	% near-deterministic embeddings (e.g., low-noise chaotic maps, whose points
	% collapse onto a thin curve), where pooledCov can be near-singular and the raw
	% separation explodes to values in the thousands despite no multi-modal
	% structure. The log keeps sep_mahal monotonic in separation while preventing
	% these degenerate cases from dominating downstream analyses:
	out.sep_mahal = log1p(sqrt(dMu / pooledCov * dMu'));
end

% ------------------------------------------------------------------------------
% valley: depth of the density valley at the sep_pp cut
% ------------------------------------------------------------------------------
out.valley = ValleyDepth(ppProj, thr(iBest));

% ------------------------------------------------------------------------------
% nongauss: negentropy of the whitened cloud (nearest-neighbor entropy estimate)
% ------------------------------------------------------------------------------
out.nongauss = Negentropy(white);

end
% ------------------------------------------------------------------------------
function [r2, thr] = OptimalTwoMeansSplits(proj)
	% Exact optimal two-cluster (1-D 2-means) split of each column of proj: the cut
	% between distinct consecutive sorted values that maximizes the between-cluster
	% sum of squares. r2 is that sum of squares as a fraction of the total (0 for a
	% constant column); thr is the cut point, midway between the values either side.
	[n, d] = size(proj);
	zs = sort(proj, 1);
	csum = cumsum(zs, 1);
	nLo = (1:n - 1)';
	meanLo = csum(1:end - 1, :) ./ nLo;
	meanHi = (csum(end, :) - csum(1:end - 1, :)) ./ (n - nLo);
	bss = nLo .* (n - nLo) / n .* (meanLo - meanHi).^2;
	bss(zs(1:end - 1, :) == zs(2:end, :)) = -Inf; % (no cut between tied values)
	[bssMax, j] = max(bss, [], 1);
	tss = sum((zs - mean(zs, 1)).^2, 1);
	r2 = zeros(d, 1);
	ok = tss > 0 & isfinite(bssMax);
	r2(ok) = bssMax(ok) ./ tss(ok);
	r2 = min(max(r2, 0), 1);
	idx = sub2ind([n, d], j, 1:d);
	thr = 0.5 * (zs(idx) + zs(idx + 1))';
end
% ------------------------------------------------------------------------------
function dirs = SplitDirections(m)
	% Fixed grid of unit directions in m dimensions (one per row; a direction and its
	% negative describe the same split, so only a hemisphere is covered)
	if m == 1
		dirs = 1;
	elseif m == 2
		theta = (0:89)' * pi / 90;
		dirs = [cos(theta), sin(theta)];
	elseif m == 3
		numDirs = 300;
		k = (0:numDirs - 1)' + 0.5;
		z = k / numDirs;
		phi = pi * (1 + sqrt(5)) * k;
		r = sqrt(1 - z.^2);
		dirs = [r .* cos(phi), r .* sin(phi), z];
	else
		% Generalized golden-ratio sequence (x^(m+1) = x + 1) mapped to the sphere
		numDirs = 300;
		phiM = 2;
		for it = 1:50
			phiM = (1 + phiM)^(1 / (m + 1));
		end
		alpha = (1 / phiM).^(1:m);
		u = mod(0.5 + (1:numDirs)' * alpha, 1);
		u = min(max(u, 1e-12), 1 - 1e-12);
		dirs = -sqrt(2) * erfcinv(2 * u); % (standard normal quantiles)
		dirs = dirs ./ sqrt(sum(dirs.^2, 2));
		dirs = dirs .* (1 - 2 * (dirs(:, 1) < 0));
	end
end
% ------------------------------------------------------------------------------
function gm = GMM2FromLabels(Y, labels, regVal)
	% Two-component full-covariance Gaussian mixture fitted by EM, started from the
	% hard assignment labels (true/false): the two groups' proportions, means and
	% regularized covariances. Each iteration is an E-step followed by an M-step
	% (covariances regularized by adding regVal to the diagonal); iterations stop
	% when the mean log-likelihood per point (evaluated at the E-step) changes by
	% less than 1e-6, or after 500 iterations. Returns a structure with fields mu,
	% Sigma and meanLogLik (the mean log-likelihood per point of the final fit), or
	% empty if a group has no more than m points or a covariance is not positive
	% definite.
	[n, m] = size(Y);
	groups = {~labels, labels};
	if min(cellfun(@nnz, groups)) <= m
		gm = []; return
	end
	mu = zeros(2, m);
	Sigma = zeros(m, m, 2);
	w = zeros(1, 2);
	for g = 1:2
		Yg = Y(groups{g}, :);
		mu(g, :) = mean(Yg, 1);
		Sigma(:, :, g) = cov(Yg) + regVal * eye(m);
		w(g) = mean(groups{g});
	end
	maxIter = 500;
	tol = 1e-6;
	lowerBound = -Inf;
	for it = 1:maxIter
		% E-step:
		[logProb, ok] = WeightedLogProb(Y, w, mu, Sigma);
		if ~ok
			gm = []; return
		end
		logProbNorm = LogSumExp(logProb);
		resp = exp(logProb - logProbNorm);
		% M-step:
		nk = sum(resp, 1) + 10 * eps;
		mu = (resp' * Y) ./ nk';
		for k = 1:2
			d = Y - mu(k, :);
			Sigma(:, :, k) = ((resp(:, k) .* d)' * d) / nk(k) + regVal * eye(m);
		end
		w = nk / sum(nk);
		% Convergence:
		previousBound = lowerBound;
		lowerBound = mean(logProbNorm);
		if abs(lowerBound - previousBound) < tol
			break
		end
	end
	[logProb, ok] = WeightedLogProb(Y, w, mu, Sigma);
	if ~ok
		gm = []; return
	end
	gm = struct('mu', mu, 'Sigma', Sigma, 'meanLogLik', mean(LogSumExp(logProb)));
end
% ------------------------------------------------------------------------------
function [logProb, ok] = WeightedLogProb(Y, w, mu, Sigma)
	% log(w_k) + log N(y | mu_k, Sigma_k) for each point (row) and component (column)
	[n, m] = size(Y);
	K = numel(w);
	logProb = zeros(n, K);
	ok = true;
	for k = 1:K
		[R, p] = chol(Sigma(:, :, k)); % (Sigma = R'*R)
		if p > 0
			ok = false; return
		end
		z = (Y - mu(k, :)) / R;
		logProb(:, k) = -0.5 * (m * log(2*pi) + sum(z.^2, 2)) - sum(log(diag(R))) + log(w(k));
	end
end
% ------------------------------------------------------------------------------
function s = LogSumExp(a)
	% log(sum(exp(a), 2)), computed stably
	aMax = max(a, [], 2);
	s = aMax + log(sum(exp(a - aMax), 2));
end
% ------------------------------------------------------------------------------
function v = ValleyDepth(proj, thr)
	% Depth of the density valley at the cut thr of the projected points proj (unit
	% variance): one minus the Gaussian-kernel density at thr relative to the smaller
	% of the densities at the two groups' lower median points, floored at 0. The
	% bandwidth is Silverman's rule of thumb for unit variance, 0.9*N^(-1/5).
	h = 0.9 * numel(proj)^(-1/5);
	kde = @(x) mean(exp(-0.5 * ((x - proj) / h).^2));
	lo = sort(proj(proj <= thr));
	hi = sort(proj(proj > thr));
	if isempty(lo) || isempty(hi)
		v = 0; return
	end
	fRef = min(kde(lo(ceil(numel(lo) / 2))), kde(hi(ceil(numel(hi) / 2))));
	v = max(0, 1 - kde(thr) / fRef);
end
% ------------------------------------------------------------------------------
function ng = Negentropy(Y)
	% Negentropy (nats) of the whitened point cloud Y (one point per row): the entropy
	% of the identity-covariance Gaussian, (m/2)*log(2*pi*e), minus the
	% Kozachenko-Leonenko nearest-neighbor entropy estimate (k = 5 neighbors), with
	% each nearest-neighbor distance floored at 0.01.
	[n, m] = size(Y);
	k = 5;
	[~, dist] = knnsearch(Y, Y, 'K', k + 1); % (each point is among its own nearest neighbors)
	r = max(dist(:, end), 0.01);
	logVolUnitBall = (m / 2) * log(pi) - gammaln(m / 2 + 1);
	H = psi(n) - psi(k) + logVolUnitBall + m * mean(log(r));
	ng = (m / 2) * log(2 * pi * exp(1)) - H;
end
