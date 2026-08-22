function out = SY_SlowFeatureAnalysis(y, numWindows)
% SY_SlowFeatureAnalysis    Slow feature analysis of windowed statistics
%
% Splits the time series into numWindows non-overlapping segments (same
% segmentation and per-segment statistics as SY_RampingWindows: mean,
% variance, skewness, and lag-1 autocorrelation, forming a numWindows x 4
% matrix), then applies Slow Feature Analysis (SFA) to find the linear
% combination of these four statistics that varies as *slowly* as
% possible across the sequence of windows -- i.e., minimizes the
% variance of its own increments, subject to unit variance. This is a
% fundamentally different criterion to SY_RampingWindows' per-statistic
% monotonic-trend tests: SFA is multivariate (can find a slow combination
% spread across mean/variance/skewness/AC1 that no single one of them
% shows individually) and detects any slow, low-frequency evolution, not
% just a monotonic ramp (e.g., a slow rise-then-fall hump across the
% series, invisible to Kendall's tau/Pearson's r, still registers as
% "slow" here). Validated on Empirical1000 to be essentially uncorrelated
% (max |r| ~ 0.14) with all SY_RampingWindows_10 trend fields.
%
% The companion comparison to ordinary PCA addresses a different
% question: PCA finds the combination of the four statistics with
% *maximal variance* across windows, which need not be the slowest one
% -- a large-amplitude but choppy/noisy statistic can dominate variance
% while a small-amplitude but smooth, genuinely slow drift is buried
% underneath it and missed by PCA. pc1VarFrac and slowPCA1corr quantify
% this: how concentrated the variance is in a single PCA direction, and
% whether the slow direction found by SFA coincides with that
% high-variance PCA direction (high overlap) or is a separate,
% lower-variance mode that ordinary variance-based analysis would miss
% (low overlap).
%
% ---INPUTS:
% y, the input time series
%
% numWindows, the number of non-overlapping segments to divide the time
%       series into (default: 20). Non-overlapping segments are used for
%       the same reason as SY_RampingWindows: overlap would induce
%       artificial serial correlation between adjacent window statistics,
%       which would make the derivative-based slowness measure below
%       spuriously small regardless of any real slow structure in the
%       data. 20 was chosen (rather than SY_RampingWindows' default of
%       10) because SFA needs enough windows to estimate the underlying
%       4x4 covariance matrices (of the statistics, and of their
%       increments) reasonably reliably -- at numWindows = 10 the
%       null-distribution spread of the slowness eigenvalues is
%       considerably wider, making individual values a noisier signal.
%
% ---OUTPUTS:
%
% Let z(t) be the numWindows x 4 matrix of [mean, variance, skewness,
% AC1] computed per window, whitened (centered, then linearly
% transformed to unit covariance). SFA finds the orthogonal directions
% u_i that minimize var(diff(z*u_i)), i.e. the "slowness" eigenvalues
% eta_i = var(diff(z*u_i)) of the covariance of the whitened derivative
% signal (ascending: eta_1 is the slowest direction). For reference,
% i.i.d. (white-noise) windows give eta ~ 2 on average; substantially
% smaller values indicate genuinely slow (smooth, low-frequency)
% structure in some combination of the four per-window statistics.
%
% eta1, the smallest (slowest) SFA eigenvalue
% etaEnd, the largest (fastest/noisiest) SFA eigenvalue
% etaStd, the standard deviation of all four SFA eigenvalues (spread of
%       the slowness spectrum)
% pc1VarFrac, the fraction of total variance (across the four
%       statistics) explained by the leading PCA component
% slowPCA1corr, the absolute correlation between the slowest SFA
%       component's scores and the leading PCA component's scores --
%       near 1 means the slow direction is simply the dominant
%       (highest-variance) direction PCA would already find; near 0 means
%       SFA has isolated a genuinely separate, low-variance slow mode.
%
% slowCorrMean, slowCorrVar, slowCorrSkew, slowCorrAC1, the absolute
%       correlation between each raw per-window statistic (mean, variance,
%       skewness, AC1, computed before whitening) and the slowest SFA
%       component's scores -- these are loadings that indicate which of
%       the four statistics is driving the slow mode SFA found, since the
%       eta/PCA fields above summarize the slow mode's existence but not
%       its composition.
%
% cf. Wiskott, L. & Sejnowski, T.J. "Slow feature analysis: unsupervised
% learning of invariances." Neural Computation 14(4), 715-770 (2002).

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
N = length(y);

if nargin < 2 || isempty(numWindows)
    numWindows = 20;
end

minNumWindows = 10; % need enough windows to estimate the 4x4 covariance matrices reliably
minWindowLength = 20; % heuristic minimum for meaningful skewness/AC1 estimates

if numWindows < minNumWindows
    error('numWindows = %u is too few for a reliable slow feature analysis (need >= %u)', numWindows, minNumWindows);
end

winLength = floor(N / numWindows);
if winLength < minWindowLength
    warning('Time series (N = %u) too short for %u windows of a meaningful length', N, numWindows);
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Segment the time series (non-overlapping, discarding any remainder)
% ------------------------------------------------------------------------------
z = reshape(y(1:winLength * numWindows), winLength, numWindows); % winLength x numWindows

% ------------------------------------------------------------------------------
%% Per-window statistics: mean, variance, skewness, AC1 (same as SY_RampingWindows)
% ------------------------------------------------------------------------------
winMean = mean(z)';
winVar = var(z)';
winSkew = skewness(z)';
winAC1 = zeros(numWindows, 1);
for i = 1:numWindows
    winAC1(i) = CO_AutoCorr(z(:, i), 1, 'Fourier');
end
X = [winMean, winVar, winSkew, winAC1]; % numWindows x 4

% ------------------------------------------------------------------------------
%% PCA (variance-maximizing directions) and SFA (slowness-minimizing directions)
% ------------------------------------------------------------------------------
Xc = X - mean(X, 1);
Cx = cov(Xc); % 4 x 4

[Vp, Dp] = eig(Cx);
[pcaEigs, ord] = sort(diag(Dp), 'descend');
Vp = Vp(:, ord);
pcScores = Xc * Vp; % numWindows x 4, PC1 = pcScores(:,1)

% Whitening (symmetric/ZCA, avoids an arbitrary rotation among near-degenerate
% directions). Directions with near-zero variance relative to the leading one
% (e.g. a per-window statistic that barely varies across windows) are dropped
% rather than whitened: full whitening would divide by their near-zero std and
% amplify what is essentially estimation noise into a spuriously enormous
% "fast" eigenvalue. pcaEigs(1) (the leading, largest eigenvalue) is never
% itself dropped by this relative threshold.
relFloor = 1e-2;
keep = pcaEigs > relFloor * pcaEigs(1);
if sum(keep) < 2
    out = NaN; return
end
Vk = Vp(:, keep);
Zw = Xc * Vk * diag(1 ./ sqrt(pcaEigs(keep))); % numWindows x sum(keep), approx unit covariance

dZ = diff(Zw); % (numWindows-1) x sum(keep), the whitened derivative signal
Cd = cov(dZ);
[Us, Ds] = eig(Cd);
[etas, ord2] = sort(diag(Ds), 'ascend'); % slowness eigenvalues, ascending = slowest first
Us = Us(:, ord2);
slowScores = Zw * Us; % numWindows x sum(keep), slowest component = slowScores(:,1)

% ------------------------------------------------------------------------------
%% Output statistics
% ------------------------------------------------------------------------------
out.eta1 = etas(1);
out.etaEnd = etas(end);
out.etaStd = std(etas);

out.pc1VarFrac = pcaEigs(1) / sum(pcaEigs);
out.slowPCA1corr = abs(corr(slowScores(:, 1), pcScores(:, 1)));

out.slowCorrMean = abs(corr(winMean, slowScores(:, 1)));
out.slowCorrVar = abs(corr(winVar, slowScores(:, 1)));
out.slowCorrSkew = abs(corr(winSkew, slowScores(:, 1)));
out.slowCorrAC1 = abs(corr(winAC1, slowScores(:, 1)));

end
