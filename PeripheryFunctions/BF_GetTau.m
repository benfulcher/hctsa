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
%           tau = 1). If the ACF never crosses 1/e, the AMI is searched up to
%           N/10 lags, taking its first local minimum, or else its first
%           crossing below -0.5*log(1 - exp(-2)) (the Gaussian AMI at a
%           correlation of 1/e).
%   (iii) 'mi-gaussian': the first local minimum of the Gaussian AMI
%           (= CO_FirstMin(y,'mi-gaussian')). The Gaussian AMI is a monotonic
%           function of |ACF|, so this is not a nonlinear timescale. This was
%           the rule used for 'mi' delays in earlier versions of hctsa.
%   (iv) 'ac': the first zero crossing of the ACF (as CO_FirstCrossing(y,'ac',0,'discrete')).
%
%---OUTPUT:
% tau, the time delay (NaN if it cannot be determined, e.g., for a constant
%       series or one whose ACF never crosses the required threshold).

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
    % min(first minimum of the Kraskov AMI, floor of the ACF 1/e crossing):
    N = length(y);
    tauAC = TauAC1e(y);
    if isnan(tauAC)
        maxLag = floor(N/10);
    else
        % (only lags up to tauAC can matter; one extra lag to detect a minimum at tauAC)
        maxLag = min(tauAC + 1, floor(N/10));
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
    amiStruct = IN_AutoMutualInfo(y, 1:maxLag, 'kraskov1', '4');
    amis = arrayfun(@(l) amiStruct.(sprintf('ami%u', l)), 1:maxLag);

    % First local minimum (a minimum at lag 1 if the AMI already increases from lag 1 to 2):
    tauMin = NaN;
    for l = 1:maxLag - 1
        if isnan(amis(l + 1))
            break
        end
        if amis(l + 1) > amis(l) && (l == 1 || amis(l - 1) > amis(l))
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
