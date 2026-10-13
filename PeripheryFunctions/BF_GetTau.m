function tau = BF_GetTau(y, tauRule)
% BF_GetTau   A time delay (in samples) set by the time series' own timescale
%
% Adaptive time delays let a feature measure structure relative to the
% series' own correlation time, rather than relative to the sampling interval
% (a fixed lag of 1 mostly measures smoothness when a process is oversampled).
% Both rules below scale linearly with the sampling rate, and are chosen to sit
% on the low side of the correlation time: iterated maps, which decorrelate
% within a single step, keep tau = 1 (so their one-step recurrence is retained).
%
%---INPUTS:
% y, the input time series (column vector)
% tauRule, how to set the delay:
%   (i) 'ac1e': the largest integer lag at which the autocorrelation is still
%           at least 1/e (the floor of the linearly interpolated first 1/e
%           crossing of the ACF), and at least 1. Decimating by this delay never
%           goes past decorrelation: the decimated series has a lag-1
%           autocorrelation of at least 1/e. The 1/e crossing is much more
%           stable than the first zero crossing ('ac'), which for a
%           monotonically decaying ACF is set mostly by estimation noise.
%   (ii) 'mi': the smaller of the first local minimum of the (Kraskov, k = 4)
%           automutual information (AMI) and the 'ac1e' delay, and at least 1.
%           The first AMI minimum is the classic Fraser-Swinney delay for
%           (nonlinear) flows; taking the smaller of the nonlinear and linear
%           decorrelation times keeps the delay low (a chaotic map has a slowly
%           decaying AMI but an ACF that decorrelates within a step, so gets
%           tau = 1). The AMI is computed on a tie-aware rank transform of the
%           series (see NOTES), and the first local minimum is taken of its
%           centered 3-point moving average (2-point at the ends), over lags
%           1, ..., ac1e + 2. If the ACF never crosses 1/e, the AMI is searched up
%           to N/10 lags, taking the first local minimum of the smoothed AMI, or
%           else the first crossing of the (unsmoothed) AMI below
%           -0.5*log(1 - exp(-2)) (the Gaussian AMI at a correlation of 1/e).
%   (iii) 'mi-gaussian': the first local minimum of the Gaussian AMI
%           (= CO_FirstMin(y,'mi-gaussian')). The Gaussian AMI is a monotonic
%           function of |ACF|, so this is not a nonlinear timescale. This was
%           the rule used for 'mi' delays in earlier versions of hctsa.
%   (iv) 'ac': the first zero crossing of the ACF (as CO_FirstCrossing(y,'ac',0,'discrete')).
%
%---OUTPUT:
% tau, the time delay (NaN if it cannot be determined, e.g., for a constant
%       series or one whose ACF never crosses the required threshold).
%
%---NOTES:
% The 'mi' delay is designed to be insensitive to rounding-level changes in the
% values and to how ties are broken. The AMI is invariant to monotone
% transformations, so it is computed on average ranks in which repeated values
% are ties: exact repeats, and also isolated clusters of near-identical values
% (a maximal run of sorted values whose consecutive gaps are all at most
% 1e-5*std(y), separated from the rest by gaps of more than 1e-4*std(y) on both
% sides). The ranks are z-scored and a fixed, reproducible jitter of
% 1e-10 x (standard normal numbers from BF_Random with seed 0) is always added
% to break the remaining ties between equal values and equal distances. The
% first minimum of a Kraskov AMI curve is often shallow and comparable to the
% estimator's own fluctuations between neighboring lags, which the 3-point
% smoothing suppresses.

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

% Many operations resolve the same delay for the same series, so cache the last
% few results (the Kraskov AMI used by 'mi' is the expensive part):
persistent cacheY cacheRule cacheTau
maxCacheEntries = 4;
if isempty(cacheY)
    cacheY = {};
    cacheRule = {};
    cacheTau = {};
end
for ci = 1:numel(cacheY)
    if isequal(tauRule, cacheRule{ci}) && isequaln(y, cacheY{ci})
        tau = cacheTau{ci};
        return
    end
end

switch tauRule
    case 'ac1e'
        tau = TauAC1e(y);
    case 'mi'
        tau = TauMI(y);
    case 'mi-gaussian'
        tau = CO_FirstMin(y, 'mi-gaussian');
    case 'ac'
        tau = CO_FirstCrossing(y, 'ac', 0, 'discrete');
    otherwise
        error('Unknown time-delay rule ''%s''', tauRule);
end

cacheY{end + 1} = y;
cacheRule{end + 1} = tauRule;
cacheTau{end + 1} = tau;
if numel(cacheY) > maxCacheEntries
    cacheY(1) = [];
    cacheRule(1) = [];
    cacheTau(1) = [];
end

end
% ------------------------------------------------------------------------------
function tau = TauAC1e(y)
    % Floor of the first 1/e crossing of the ACF (at least 1):
    acf = CO_AutoCorr(y, [], 'Fourier'); % lags 0, 1, ..., N-1
    threshold = 1/exp(1);
    if any(isnan(acf)) || ~any(acf < threshold)
        % Degenerate series, or the ACF never decays to 1/e
        tau = NaN;
        return
    end
    [~, pointOfCrossing] = BF_PointOfCrossing(acf, threshold);
    tau = max(1, floor(pointOfCrossing - 1)); % (index -> lag)
end
% ------------------------------------------------------------------------------
function tau = TauMI(y)
    % min(first minimum of the smoothed Kraskov AMI of tie-aware ranks, floor of the
    % ACF 1/e crossing):
    N = length(y);
    tauAC = TauAC1e(y);
    if isnan(tauAC)
        maxLag = floor(N/10);
    else
        % (only lags up to tauAC can matter; two extra lags so the smoothed curve can
        % show a minimum at tauAC)
        maxLag = min(tauAC + 2, floor(N/10));
    end
    if tauAC == 1
        % Can't go below 1
        tau = 1;
        return
    end
    if maxLag < 2
        % Time series too short to locate an AMI minimum
        tau = tauAC;
        return
    end

    % Tie-aware ranks, z-scored, plus a fixed reproducible jitter:
    r = TieAwareRanks(y(:));
    if ~(std(r) > 0)
        tau = NaN; % (all values tied)
        return
    end
    r = (r - mean(r)) / std(r);
    r = r + 1e-10 * BF_Random(N, 0, 'normal');

    amiStruct = IN_AutoMutualInfo(r, 1:maxLag, 'kraskov1', '4');
    amis = arrayfun(@(l) amiStruct.(sprintf('ami%u', l)), 1:maxLag);
    numValid = find(isnan(amis), 1, 'first') - 1; % (NaNs, for too few samples, only trail)
    if ~isempty(numValid)
        amis = amis(1:numValid);
    end
    if numel(amis) < 2
        tau = tauAC;
        return
    end

    % Centered 3-point moving average (2-point at the ends):
    amisSmooth = amis;
    amisSmooth(2:end - 1) = (amis(1:end - 2) + amis(2:end - 1) + amis(3:end)) / 3;
    amisSmooth(1) = (amis(1) + amis(2)) / 2;
    amisSmooth(end) = (amis(end - 1) + amis(end)) / 2;

    % First local minimum of the smoothed curve (at lag 1 if it already increases
    % from lag 1 to 2):
    tauMin = NaN;
    for l = 1:numel(amisSmooth) - 1
        if amisSmooth(l + 1) > amisSmooth(l) && (l == 1 || amisSmooth(l - 1) > amisSmooth(l))
            tauMin = l;
            break
        end
    end

    if ~isnan(tauMin)
        tau = min(tauMin, tauAC); % (min ignores a NaN tauAC)
    elseif ~isnan(tauAC)
        tau = tauAC;
    else
        % No 1/e crossing of the ACF and no AMI minimum: the first drop of the AMI
        % below the Gaussian AMI at a correlation of 1/e
        miThreshold = -0.5*log(1 - exp(-2));
        tau = find(amis < miThreshold, 1, 'first');
        if isempty(tau)
            tau = NaN;
        end
    end
end
% ------------------------------------------------------------------------------
function r = TieAwareRanks(y)
    % Average ranks in which exact repeats, and isolated clusters of near-identical
    % values, are ties. A maximal run of sorted values with consecutive gaps of at
    % most tol*std(y) is one tied value if the gaps before and after the run both
    % exceed 10*tol*std(y); otherwise only exactly equal values within it are tied.
    tol = 1e-5;
    iso = 10;
    n = numel(y);
    [ys, idx] = sort(y);
    g = diff(ys);
    t = tol * std(y);
    newRun = [true; g > t];
    runID = cumsum(newRun);
    runCount = accumarray(runID, 1);
    runLast = cumsum(runCount);
    runFirst = runLast - runCount + 1;
    gPad = [Inf; g; Inf]; % gPad(k): gap before sorted point k; gPad(k + 1): gap after it
    isolated = (gPad(runFirst) > iso * t) & (gPad(runLast + 1) > iso * t);
    newGroup = [true; g > 0];
    inIsolated = isolated(runID);
    newGroup(inIsolated) = newRun(inIsolated);
    groupID = cumsum(newGroup);
    groupCount = accumarray(groupID, 1);
    groupStart = cumsum(groupCount) - groupCount;
    midRank = groupStart + (groupCount + 1) / 2;
    r = zeros(n, 1);
    r(idx) = midRank(groupID);
end
