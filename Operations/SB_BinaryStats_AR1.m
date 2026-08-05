function out = SB_BinaryStats_AR1(y, binaryMethod)
% SB_BinaryStats_AR1    Binary run-length statistics normalized against an AR(1) null
%
% Binarizes the time series (as in SB_BinaryStats) and compares the
% resulting run-length statistics to their analytic expectation under a
% Gaussian AR(1) null process with the same lag-1 autocorrelation as y.
%
% ---INPUTS:
% y, the input time series
%
% binaryMethod, the symbolization rule -- 'mean' (above/below the mean) or
%               'diff' (sign of incremental differences). Unlike
%               SB_BinaryStats, 'iqr' is not supported here: the theory
%               below is specific to a sign/half-line threshold at the
%               mean of a (possibly transformed) Gaussian series, and does
%               not apply to interval/band membership. Default: 'mean'.
%
% ---THEORY:
% For a stationary Gaussian process u (u = y for 'mean'; u = diff(y) for
% 'diff', i.e., whichever series BF_Binarize actually thresholds at zero),
% the probability that consecutive samples lie on the same side of the
% mean follows the classical arcsine law for bivariate normal sign
% agreement:
%
%   p = P(sign(u_t - mean(u)) == sign(u_{t+1} - mean(u))) = 1/2 + asin(rho)/pi
%
% where rho is the lag-1 autocorrelation of u. This is exact for any
% bivariate normal pair with correlation rho. Treating the resulting
% binary sign sequence as an (approximate) two-state Markov chain with
% persistence probability p gives geometrically distributed run lengths,
% hence:
%
%   E[run length]                     = 1/(1-p)
%   E[proportion of series in 1-runs] = (1-p)/2  (by symmetry about the mean)
%
% Because these depend only on the *marginal* (pairwise) transition
% probability, they hold up well empirically for real Gaussian AR(1)
% processes even at high persistence:
%   - 'mean': verified by simulation across phi = 0-0.95 (empirical-to-
%     expected ratio within ~1% of 1 throughout).
%   - 'diff': verified similarly across phi = -0.85-0.95 (ratio within
%     ~0.2% of 1 throughout -- even tighter than 'mean', since
%     differencing keeps the diff series' own lag-1 autocorrelation away
%     from +-1 even as phi -> 1).
%
% The same two-state-Markov assumption also yields closed forms for the
% variance of run lengths and for low-order point probabilities (e.g.,
% P(L=2)-P(L=1)), but these do NOT hold up empirically: simulation (for
% 'mean') shows the std-ratio drifting from ~1.0 at phi=0 to ~1.7 at
% phi=0.95, because the true sign sequence of a Gaussian AR(1) is a
% hidden- rather than first-order Markov chain (a sign flip's probability
% depends on the current value's magnitude, not just its sign, giving
% real run lengths a heavier tail than the geometric model predicts). The
% mean is robust to this mismatch; the variance and point probabilities
% are not. For the same reason (and because it's an extreme-value
% statistic besides), no closed form is attempted for the *longest* run
% either. So this operation is deliberately narrow: only meanstretch0/1
% and pstretch1 -- the statistics with theory that actually holds up --
% get an AR(1)-normalized counterpart. (For the fuller set of empirical
% run-length statistics, including longstretch0/1 and stdstretch0/1, see
% SB_BinaryStats.)
%
% ---OUTPUTS: the AR(1)-implied persistence probability and analytic
% expectations, and the empirical-to-expected ratios for meanstretch0/1 and
% pstretch1. The raw empirical meanstretch0/1 and pstretch1 values
% themselves are not separately registered as features here since they are
% identical to SB_BinaryStats_mean's/SB_BinaryStats_diff's fields of the
% same name (both use the same binarization and run-length extraction).

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
if nargin < 2 || isempty(binaryMethod)
    binaryMethod = 'mean';
end
if ~ismember(binaryMethod, {'mean', 'diff'})
    error('SB_BinaryStats_AR1 supports binaryMethod ''mean'' or ''diff'' only (not ''%s'')', binaryMethod);
end

% ------------------------------------------------------------------------------
%% Identify the series that is actually being sign-thresholded at its mean
% ------------------------------------------------------------------------------
switch binaryMethod
    case 'mean'
        u = y;
    case 'diff'
        u = diff(y);
end

% ------------------------------------------------------------------------------
%% Binarize (the case for which the arcsine law below is exact)
% ------------------------------------------------------------------------------
yBin = BF_Binarize(y, binaryMethod);
N = length(yBin); % note: N = length(y)-1 for 'diff'

% ------------------------------------------------------------------------------
%% Empirical run-length statistics (cf. SB_BinaryStats)
% ------------------------------------------------------------------------------
difffy = diff(find([1; yBin; 1]));
stretch0 = difffy(difffy ~= 1) - 1;

difffy = diff(find([0; yBin; 0] == 0));
stretch1 = difffy(difffy ~= 1) - 1;

out.pstretch1 = length(stretch1) / N;
out.meanstretch0 = mean(stretch0); % empty stretch0 (all 1s) gives mean([]) = NaN
out.meanstretch1 = mean(stretch1); % empty stretch1 (all 0s) gives mean([]) = NaN

% ------------------------------------------------------------------------------
%% AR(1)-null persistence probability, via the arcsine law
% ------------------------------------------------------------------------------
rho = CO_AutoCorr(u, 1, 'Fourier');
rho = max(min(rho, 1), -1); % guard against tiny numerical overshoot outside [-1,1]
p = 0.5 + asin(rho) / pi;
out.ar1_p = p; % the AR(1)-implied persistence probability itself

if p >= 1 - 1e-8 % degenerate limit (rho -> 1): expected run length diverges
    out.meanstretch_ar1exp = Inf;
    out.meanstretch0_ar1rat = NaN; out.meanstretch1_ar1rat = NaN;
    out.pstretch1_ar1exp = 0;
    out.pstretch1_ar1rat = NaN;
    return
end

expMeanStretch = 1 / (1 - p);
expPstretch1 = (1 - p) / 2;

% Note: the null expectation for meanstretch0 and meanstretch1 is the same
% value (symmetric about the mean by construction), so only one expected-
% value field is needed -- registering both would be a duplicate.
out.meanstretch_ar1exp = expMeanStretch;
out.meanstretch0_ar1rat = out.meanstretch0 / expMeanStretch;
out.meanstretch1_ar1rat = out.meanstretch1 / expMeanStretch;

out.pstretch1_ar1exp = expPstretch1;
out.pstretch1_ar1rat = out.pstretch1 / expPstretch1;

end
