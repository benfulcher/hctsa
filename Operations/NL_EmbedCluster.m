function out = NL_EmbedCluster(y, tau, m, kMax, maxN)
% NL_EmbedCluster   Gaussian-mixture clustering structure in a time-delay embedding space
%
% Reconstructs the time series as a time-delay embedding (as in NL_EmbedPCA,
% NL_EmbedKernelPCA) and fits Gaussian mixture models with a small grid of
% component counts to the resulting point cloud. A dynamical process whose
% trajectory visits distinct regions of phase space (e.g., alternating
% between two attractor states, or a system with intermittent bursts) leaves
% a multi-modal point cloud in the embedding; a process with a single smooth
% (e.g., unimodal-stochastic or single-loop periodic) attractor does not.
% This is a distinct signal from marginal-distribution multi-modality (e.g.
% DN_ kurtosis/bimodality stats on the raw values), since two states can
% overlap entirely in amplitude yet still separate cleanly once lagged
% coordinates are added, and from regime-switching detected by MF_hmm_Fit /
% MF_hmm_CompareNStates, which cluster points in raw-amplitude (not
% lagged/embedded) space.
%
% Rather than reporting only the BIC-optimal number of components (a
% discrete, model-selection-driven output that can be noisy/discontinuous
% across similar time series), the main outputs are continuous separation
% statistics from a *fixed* 2-component fit, alongside the (secondary)
% BIC-optimal k for reference.
%
% ---INPUTS:
% y, the input time series
%
% tau, the time-delay, can be an integer or 'ac', or 'mi' for first
%               zero-crossing of the autocorrelation function or first minimum
%               of the automutual information, respectively
%
% m, the embedding dimension (default: 2)
%
% kMax, the maximum number of Gaussian mixture components to consider when
%       searching for the BIC-optimal component count (default: 4)
%
% maxN, the maximum number of embedded points used to fit the mixture models.
%       Longer embeddings are reduced to their first maxN points (default:
%       2000; a warning is issued whenever this cropping actually happens,
%       since separation estimates are still noisy at this size and keep
%       sharpening with more -- this is a memory/time cap, not a convergence
%       point). Can be set to 'full' to disable.
%
% ---OUTPUTS:
% bestK, the BIC-optimal number of mixture components over 1:kMax
% dBIC, the relative BIC improvement of the best fit over a single
%       (unimodal) Gaussian fit; 0 when bestK == 1 (no clustering evidence)
% sep_mahal, log1p-compressed Mahalanobis separation between the two
%       component means of a fixed 2-component fit, using their pooled
%       covariance (log-compressed to tame the heavy tail from
%       near-singular covariance on near-deterministic embeddings; see NOTES)
% sep_conf, mean posterior cluster-assignment confidence (mean of the larger
%       of each point's two posterior probabilities) under the 2-component
%       fit; ranges [0.5, 1], with 0.5 indicating fully ambiguous assignment
% sep_silh, mean silhouette value of the hard (posterior-argmax) 2-cluster
%       assignment
% sep_weightbalance, ratio of the smaller to the larger mixture weight under
%       the 2-component fit; 1 for balanced clusters, ->0 as one component
%       comes to dominate (degenerating towards a unimodal fit)
%
% ---NOTES:
% The sep_* outputs are always computed from a fixed 2-component fit, so they
% are not gated by whether that fit is actually favoured by BIC over a single
% Gaussian: even a genuinely unimodal-but-elongated point cloud (e.g. AR(1)
% noise) gets split into two "confident" halves with non-trivial sep_conf/
% sep_mahal. Likewise, a curved-but-unimodal manifold (e.g. the ring traced
% out by a periodic signal in a 2-d embedding) is poorly fit by any single
% elliptical Gaussian and so also drives bestK/dBIC up, despite having no
% distinct dynamical states in the sense this feature was designed to catch.
% dBIC == 0 (bestK == 1) is a clean "no mixture structure at all" signal, but
% dBIC > 0 does not by itself distinguish true multi-modality from curvature;
% cross-reference against a periodicity/shape feature (e.g. SP_Summaries,
% CO_Embed2_Shapes) when that distinction matters. Confirmed empirically on a
% two-regime switching process (bestK=3-4, dBIC>1, sep_mahal~2.7, sep_silh~0.98)
% vs. AR(1) noise and white noise (bestK=1, dBIC=0, sep_mahal~0.1-0.7) vs. a
% clean sine wave (bestK=4, dBIC~0.4, sep_mahal~1.3 -- curvature, not
% multi-modality).
%
% Redundancy-checked (r>=0.9 threshold) against MF_hmm_Fit,
% MF_hmm_CompareNStates, CO_Embed2_Shapes and NL_EmbedKernelPCA on Bonn EEG
% (500 series, max|r|=0.79) and Empirical1000 (1000 series). On Empirical1000,
% the pre-log1p sep_mahal correlated r=0.98 (Pearson) with
% MF_hmm_CompareNStates.chLLtrain -- but this was entirely driven by two
% chaotic-map series (logistic, Ricker) whose near-noiseless embeddings
% collapse onto a thin curve, producing near-singular pooled covariance and
% raw Mahalanobis separations >2000; Spearman r was only 0.50, and trimming
% the top 1% by |chLLtrain| dropped Pearson r to 0.21. The log1p compression
% above (motivated independently, by the same degenerate cases) resolves this
% apparent redundancy along with the numerical issue.

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

if nargin < 4 || isempty(kMax)
	kMax = 4; % search up to 4 mixture components
end

if nargin < 5 || isempty(maxN)
	maxN = 2000; % caps the cost of repeated EM fits
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

% Need enough points, relative to m and kMax, for a well-posed full-covariance
% fit at the largest component count considered:
Nembed = size(y_embed, 1);
if Nembed < 20 * m * kMax
	fprintf(1, 'Not enough embedded points (%u) for a stable %u-dimensional, up-to-%u-component mixture fit\n', ...
			Nembed, m, kMax);
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
	y_gmm = y_embed;
else
	if Nembed > maxN
		warning(['Cropping to the first %u of %u embedded points for mixture fitting ' ...
				 '(memory/time cap, not a convergence point -- raise maxN or use ''full'' ' ...
				 'for a more precise, more expensive estimate)'], maxN, Nembed);
		y_gmm = y_embed(1:maxN, :);
	else
		y_gmm = y_embed;
	end
end
N = size(y_gmm, 1);

% Regularize covariance estimates proportionally to the data's own scale
% (y is z-scored, but the embedded coordinates can still have any variance):
regVal = 1e-6 * mean(var(y_gmm));
gmOptions = statset('MaxIter', 500, 'Display', 'off');

% ------------------------------------------------------------------------------
% Fit a grid of mixture models, k = 1, ..., kMax, and record their BIC
% ------------------------------------------------------------------------------
BIC = NaN(kMax, 1);
gmModels = cell(kMax, 1);

try
	gmModels{1} = fitgmdist(y_gmm, 1);
	BIC(1) = gmModels{1}.BIC;
catch
	% A single-component Gaussian fit failing (degenerate covariance) means
	% no useful mixture structure can be assessed either:
	out = NaN; return
end

% Some replicates legitimately fail to converge within MaxIter (the best of
% the remaining replicates is still used); this is expected often enough
% across a large batch of time series that the warning is silenced here:
warningState = warning('off', 'stats:gmdistribution:FailedToConvergeReps');
for k = 2:kMax
	try
		gmModels{k} = fitgmdist(y_gmm, k, 'CovarianceType', 'full', ...
								'RegularizationValue', regVal, 'Replicates', 3, ...
								'Options', gmOptions, 'Start', 'plus');
		BIC(k) = gmModels{k}.BIC;
	catch
		% This k failed to converge to a valid fit (e.g., an empty or
		% near-singular component); leave BIC(k) as NaN and move on
	end
end
warning(warningState);

[~, bestK] = min(BIC);
out.bestK = bestK;
if bestK == 1
	out.dBIC = 0;
else
	out.dBIC = (BIC(1) - BIC(bestK)) / abs(BIC(1));
end

% ------------------------------------------------------------------------------
% Continuous separation statistics from a fixed 2-component fit
% ------------------------------------------------------------------------------
if kMax < 2 || isempty(gmModels{2})
	out.sep_mahal = NaN;
	out.sep_conf = NaN;
	out.sep_silh = NaN;
	out.sep_weightbalance = NaN;
	return
end

gm2 = gmModels{2};
post = posterior(gm2, y_gmm);
[maxPost, clustIdx] = max(post, [], 2);

out.sep_conf = mean(maxPost);
out.sep_weightbalance = min(gm2.ComponentProportion) / max(gm2.ComponentProportion);

pooledCov = (gm2.Sigma(:, :, 1) + gm2.Sigma(:, :, 2)) / 2;
dMu = gm2.mu(1, :) - gm2.mu(2, :);
mahal_raw = sqrt(dMu / pooledCov * dMu');
% log1p-compressed: raw Mahalanobis separation is numerically unstable for
% near-deterministic embeddings (e.g. low-noise chaotic maps, whose points
% collapse onto a thin curve), where pooledCov can be near-singular and
% mahal_raw explodes to values in the thousands despite no real multi-modal
% structure (empirically confirmed on the logistic/Ricker map series in the
% Empirical1000 corpus: mahal_raw > 2000 vs. a typical well-separated
% two-regime process giving mahal_raw ~ 10-15). The log keeps sep_mahal
% monotonic in separation while preventing these degenerate cases from
% dominating any downstream (e.g. z-scored, correlation-based) analysis:
out.sep_mahal = log1p(mahal_raw);

if numel(unique(clustIdx)) < 2
	% All points collapsed onto one component under hard assignment: no
	% separation to measure
	out.sep_silh = 0;
else
	s = silhouette(y_gmm, clustIdx);
	out.sep_silh = mean(s);
end

end
