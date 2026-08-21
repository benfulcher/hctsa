function out = NL_EVTLocalDim(y, tau, m, q, theilerWin, nPoles, mOrder, maxN, randomSeed)
% NL_EVTLocalDim    Extreme value theory local dimension and persistence of the attractor.
%
% Time-delay embeds the series and, for a sample of reference points
% ("poles") on the reconstructed orbit, treats close returns of the orbit
% to each pole as extreme events: g_i = -log(||Y_i - pole||) is large
% exactly when the orbit passes close to the pole. Extreme value theory
% applied to this observable gives two local quantities per pole
% (Freitas, Freitas & Todd 2010; Lucarini et al., "Extremes and
% Recurrence in Dynamical Systems", 2016; Faranda, Messori & Yiou 2017):
%
%   - a local dimension d(pole): under the Freitas-Freitas-Todd theorem,
%     the Gumbel-law scale parameter of the extreme value law for g_i
%     equals the local dimension of the attractor at that pole exactly
%     (Faranda, Freitas, Guiraud & Vaienti 2014, Props 1-2). In practice
%     this is estimated by a peaks-over-threshold fit: take the
%     exceedances of g_i above a high quantile q, and set
%     d(pole) = 1/mean(exceedances) (the reciprocal of the exponential
%     MLE scale, i.e. the GPD fit with shape fixed at its ansatz-implied
%     value of 0 -- appropriate here because g_i = -log(distance) is
%     unbounded above, putting it in the Gumbel/exponential-tail domain).
%   - a persistence theta(pole) (the "extremal index" of g_i at that
%     pole): whether close returns to the pole arrive as isolated events
%     (theta near 1) or cluster into runs where the orbit lingers nearby
%     (theta << 1, i.e. long average residence time near that point of
%     phase space -- 1/theta is the average cluster/sojourn size).
%     Estimated with the O'Brien order-m estimator (Caby, Faranda,
%     Vaienti & Yiou, J. Stat. Phys. 2019, Eqs 19/21), which that paper's
%     own comparison against the alternative Süveges (2007) likelihood
%     estimator found more reliable for this observable, particularly
%     near near-periodic (sticky) poles.
%
% This differs mechanistically from hctsa's existing attractor-dimension
% operations (NL_Dimensions' box-counting/correlation-sum dimension;
% NL_FractalDimensions' nearest-neighbour-moment generalized dimensions):
% those pool all pairwise distances/neighbour ranks into one global
% scaling exponent, whereas this estimates a genuinely *local* dimension
% and persistence separately at each of several poles and reports how
% they're distributed (and covary) across the attractor -- capturing
% multifractal-style local heterogeneity that a single global exponent
% cannot.
%
% ---INPUTS:
% y, the input time series (assumed z-scored)
%
% tau, embedding time delay fed to BF_Embed (default: 'ac', first
%      zero-crossing of the autocorrelation function)
%
% m, embedding dimension fed to BF_Embed (default: 3)
%
% q, the quantile level defining "extreme" close returns: exceedances of
%    g_i = -log(distance) above its q-quantile are treated as events
%    (default: 0.98, i.e. the closest 2% of returns to each pole)
%
% theilerWin, Theiler window excluding temporally-correlated neighbours of
%             each pole from being treated as (trivially close) returns
%             (a proportion of the embedded length if in (0,1); default: 0.01)
%
% nPoles, number of reference points (poles) to sample from the embedded
%         orbit (cost is O(nPoles*Nemb); default: 200)
%
% mOrder, the order m of the O'Brien persistence estimator -- how many
%         steps ahead to check for a further exceedance before counting a
%         given exceedance as "isolated" (default: 5, following Caby et al.)
%
% maxN, maximum number of samples to consider. Defaults to 'full' (no
%       cropping) -- cost is O(nPoles*Nemb) and measured flat at ~0.3-0.4s
%       from N=10000 to N=20000, cheap enough not to crop by default; a
%       warning is given above N = 50000 if it is ever set explicitly.
%       Set to a number to crop longer series to their first maxN points.
%
% randomSeed, whether (and how) to reset the random seed, using
%             BF_ResetSeed, before sampling the poles (default: 'default')
%
% ---OUTPUTS: mean and standard deviation of the local dimension and
% persistence across poles, their correlation across poles, and the
% proportion of poles that yielded a valid estimate (too few exceedances
% otherwise -- a diagnostic of whether q/nPoles/series length were
% adequate, not a property of the dynamics itself).

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

if nargin < 2 || isempty(tau)
    tau = 'ac';
end
if nargin < 3 || isempty(m)
    m = 3;
end
if nargin < 4 || isempty(q)
    q = 0.98;
end
if nargin < 5 || isempty(theilerWin)
    theilerWin = 0.01;
end
if nargin < 6 || isempty(nPoles)
    nPoles = 200;
end
if nargin < 7 || isempty(mOrder)
    mOrder = 5;
end
if nargin < 8 || isempty(maxN)
    maxN = 'full';
end
if nargin < 9 || isempty(randomSeed)
    randomSeed = 'default';
end

if ischar(maxN) && strcmp(maxN, 'full')
    slowThreshold = 50000;
    if N > slowThreshold
        warning('Time series (%u samples) exceeds %u with maxN=''full''; computation may be slow', N, slowThreshold);
    end
elseif N > maxN
    warning('Time series (%u samples) exceeds maxN = %u; analyzing the first %u samples', N, maxN, maxN);
    y = y(1:maxN);
    N = length(y);
end

% ------------------------------------------------------------------------------
%% Embed the signal
% ------------------------------------------------------------------------------
Y = BF_Embed(y, tau, m, false);
if isscalar(Y) && isnan(Y) % embedding failed
    out = NaN; return
end
Nemb = size(Y, 1);

if (theilerWin > 0) && (theilerWin < 1) % specify a proportion
    theilerWin = round(theilerWin * Nemb);
end

minExceed = 15; % minimum exceedances required for a pole's estimate to be trusted
minNemb = ceil(minExceed / (1 - q)) + 2 * theilerWin + mOrder + 1;
if Nemb < max(minNemb, 100)
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Sample poles (reference points) from the embedded orbit
% ------------------------------------------------------------------------------
BF_ResetSeed(randomSeed);
nPoles = min(nPoles, Nemb);
poleIdx = randperm(Nemb, nPoles);

localDim = nan(nPoles, 1);
theta = nan(nPoles, 1);

for p = 1:nPoles
    j = poleIdx(p);

    % Exclude the Theiler window around this pole, keeping the rest of the
    % orbit in its original chronological order:
    keepIdx = setdiff(1:Nemb, max(1,j-theilerWin):min(Nemb,j+theilerWin));

    dist_j = sqrt(sum((Y(keepIdx, :) - Y(j, :)).^2, 2));
    isValid = dist_j > 0; % excludes exact duplicate embedded points
    if sum(isValid) < minExceed / (1 - q)
        continue
    end
    g = -log(dist_j(isValid)); % chronological order preserved (keepIdx is sorted)

    u = quantile(g, q);
    E = g > u; % exceedance indicator, in chronological order
    Nu = sum(E);
    if Nu < minExceed
        continue
    end

    % Local dimension: reciprocal of the mean exceedance (exponential/GPD-
    % with-shape-0 scale MLE), per the Freitas-Freitas-Todd theorem:
    localDim(p) = 1 / mean(g(E) - u);

    % Persistence: O'Brien order-mOrder estimator (Caby et al. 2019, Eq 19),
    % vectorized via a cumulative sum of the exceedance indicator.
    Ntot = length(E);
    if Ntot > mOrder + 1
        cumE = [0; cumsum(E)];
        validI = (1:(Ntot - mOrder))';
        futureSum = cumE(validI + mOrder + 1) - cumE(validI + 1);
        isIsolated = E(validI) & (futureSum == 0);
        numer = sum(isIsolated) / (Ntot - mOrder);
        denom = Nu / Ntot;
        theta(p) = min(numer / denom, 1); % clip finite-sample overshoot above 1
    end
end

% ------------------------------------------------------------------------------
%% Outputs
% ------------------------------------------------------------------------------
validPoles = ~isnan(localDim) & ~isnan(theta);
out.propValidPoles = mean(~isnan(localDim));

out.meanLocalDim = mean(localDim, 'omitnan');
out.stdLocalDim = std(localDim, 'omitnan');
out.meanTheta = mean(theta, 'omitnan');
out.stdTheta = std(theta, 'omitnan');

if sum(validPoles) >= 10 && std(localDim(validPoles)) > 0 && std(theta(validPoles)) > 0
    out.corrDimTheta = corr(localDim(validPoles), theta(validPoles));
else
    out.corrDimTheta = NaN;
end

end
