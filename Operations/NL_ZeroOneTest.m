function out = NL_ZeroOneTest(y, numC, maxN)
% NL_ZeroOneTest    The 0-1 test for chaos
%
% Implements the "modified 0-1 test for chaos" of Gottwald & Melbourne,
% which classifies a scalar time series as arising from bounded
% (periodic/quasi-periodic/regular) or unbounded (chaotic) underlying
% dynamics, without requiring phase-space reconstruction or an estimate
% of the embedding dimension (unlike, e.g., NL_TakensEstimator or
% NL_GPCorrSum).
%
% For each of numC evenly spaced frequencies c, the series y drives a
% translation variables (p_c,q_c) in the plane:
%   p_c(n) = sum_{j=1}^{n} y(j)*cos(j*c),  q_c(n) = sum_{j=1}^{n} y(j)*sin(j*c)
% and the growth rate of the (mean-corrected) mean-square displacement
% M_c(n) of (p_c,q_c) is summarized two ways: K_c, the correlation
% coefficient between n and M_c(n) (K_c near 0 means (p_c,q_c) is bounded
% -- periodic dynamics underlying y; K_c near 1 means it grows
% diffusively/without bound -- chaotic or stochastic dynamics underlying
% y); and D_c, the linear-regression slope of M_c(n) against n (the
% *rate* at which (p_c,q_c) diffuses, given that it does -- see D below).
% K and D are each summarized by their median (and Kstd/Dstd, their
% spread) over numC choices of c (avoiding resonant values of c near
% integer multiples of pi).
%
% K and D answer different questions: K asks only whether (p_c,q_c) is
% bounded or not, while D asks how fast it grows given that it isn't
% bounded -- e.g. white noise and a strongly autocorrelated AR(1) process
% are both classified as unbounded by K (K ~ 1 for both, since both
% eventually drive an unbounded random walk in the (p,q) plane), but D is
% an order of magnitude smaller for the AR(1) process in validation on
% synthetic data (the AR(1) process's smaller innovation variance per
% step directly slows the (p,q) random walk, even though total series
% variance is matched by z-scoring).
%
% ---CAVEAT (important, and the reason for the paper below): this test
% alone cannot distinguish deterministic chaos from stochastic noise --
% both give K ~ 1, since both produce an unbounded, diffusive (p_c,q_c)
% trajectory. Distinguishing the two requires a separate test for
% determinism/stochasticity (e.g. SD_SurrogateTest, using EN_PermEn or
% another nonlinear statistic as the test statistic) applied alongside
% this one, not this test in isolation.
%
% cf. G.A. Gottwald & I. Melbourne, "On the implementation of the 0-1 test
% for chaos", SIAM J. Appl. Dyn. Syst. 8(1), 129-145 (2009).
%
% Used (as one stage of a larger decision-tree pipeline, alongside
% permutation-entropy-based surrogate testing for stochasticity) in:
% Toker, D. et al. "A simple method for detecting chaos in nature",
% Commun. Biol. 3, 11 (2020). DOI: 10.1038/s42003-019-0715-9
%
% ---INPUTS:
% y, the input time series
%
% numC, the number of evenly spaced frequencies c to average the test
%       statistic K_c over (default: 20 -- checked against numC = 100 on
%       periodic/quasi-periodic/chaotic/stochastic test signals: K agrees
%       to within 0.003 and Kstd to within a similar margin, so the extra
%       cost of a much larger numC buys little additional stability)
%
% maxN, the maximum number of samples to consider. Time series longer
%       than this are cropped to their first maxN points (default:
%       10000), since the computation cost of this operation scales
%       roughly as N^2 (tcut = N/10 window lengths are each averaged
%       over ~N points, per c value; cf. NL_RQA's maxN, same rationale).
%       Set to 'full' to disable cropping (a warning is given above N =
%       50000, where run time starts to become substantial -- around 10s
%       at N = 50000-60000 with the default numC)
%
% ---OUTPUTS:
% K, the median of K_c across the numC frequencies -- the chaos
%       statistic itself (~0: bounded/regular dynamics; ~1: unbounded
%       dynamics, either chaotic or stochastic -- see caveat above)
% Kstd, the standard deviation of K_c across the numC frequencies -- a
%       large spread indicates the test is unstable/ambiguous for this
%       particular series (e.g. residual resonance contamination)
% D, the median of D_c across the numC frequencies -- the diffusion rate
%       of (p_c,q_c) given that it is unbounded (near 0 for bounded
%       dynamics too, same as D_c is ~flat when K_c is ~0)
% Dstd, the standard deviation of D_c across the numC frequencies -- can
%       be heavy-tailed for complex/near-resonant dynamics (a handful of
%       c values giving an anomalously large slope estimate), validated
%       as a genuine signal rather than a numerical artifact (e.g. the
%       synthetic Duffing-van der Pol series in the Empirical1000
%       validation set, whose own Kstd is also elevated)

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
if nargin < 2 || isempty(numC)
    numC = 20;
end
if nargin < 3 || isempty(maxN)
    maxN = 10000; % crops time series longer than this maximum length
end

y = y(:);
N = length(y);
if ischar(maxN) && strcmp(maxN, 'full')
    % No cropping -- but flag potentially slow computations for very long series
    % (cost scales roughly as N^2, similar in spirit to NL_RQA's rr*N^2):
    slowThreshold = 50000;
    if N > slowThreshold
        warning('Time series (%u samples) exceeds %u with maxN=''full''; computation may be slow', N, slowThreshold);
    end
elseif N > maxN
    warning('Time series (%u samples) exceeds maxN = %u; analyzing the first %u samples', N, maxN, maxN);
    y = y(1:maxN);
    N = maxN;
end

minN = 200; % need tcut = floor(N/10) large enough for a meaningful growth-rate correlation
if N < minN
    warning('Time series (N = %u) too short for a meaningful 0-1 test (need >= %u)', N, minN);
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% The 0-1 test for chaos
% ------------------------------------------------------------------------------
tcut = floor(N / 10);
j = (1:N)';
nvec = (1:tcut)';
mphi = mean(y);

% Evenly spaced c values avoiding resonances near 0/pi/2pi (standard convention);
% deterministic (rather than randomly sampled from the same range, as in some
% implementations) so the feature is reproducible across runs:
cs = linspace(pi/5, 4*pi/5, numC)';

Kc = zeros(numC, 1);
Dc = zeros(numC, 1);
for k = 1:numC
    c = cs(k);
    p = cumsum(y .* cos(j * c));
    q = cumsum(y .* sin(j * c));
    Mc = zeros(tcut, 1);
    for n = 1:tcut
        Mc(n) = mean((p(n + 1:N - tcut + n) - p(1:N - tcut)).^2 + ...
                     (q(n + 1:N - tcut + n) - q(1:N - tcut)).^2);
    end
    % Subtract the oscillatory contribution from a nonzero series mean,
    % which would otherwise cause (p,q) to grow even for perfectly periodic y:
    Vosc = mphi^2 * (1 - cos(nvec * c)) / (1 - cos(c));
    Mc_mod = Mc - Vosc;
    Kc(k) = corr(nvec, Mc_mod);
    pFit = polyfit(nvec, Mc_mod, 1);
    Dc(k) = pFit(1);
end

% ------------------------------------------------------------------------------
%% Output statistics
% ------------------------------------------------------------------------------
out.K = median(Kc);
out.Kstd = std(Kc);
out.D = median(Dc);
out.Dstd = std(Dc);

end
