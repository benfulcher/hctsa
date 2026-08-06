function out = SP_Bicoherence(y, segLength, maxN, numSurr)
% SP_Bicoherence   Squared bicoherence of a time series
%
% Estimates the bicoherence, a normalized bispectrum, using the standard
% segment-averaging (indirect FFT) method: the series is split into
% overlapping segments (50% overlap, as many as fit), each segment is
% demeaned, Hamming-windowed, and Fourier transformed, and the resulting
% per-segment Fourier coefficients are combined into a triple-product
% estimate of the bispectrum B(f1,f2) = <X(f1) X(f2) X*(f1+f2)>, averaged
% over segments and normalized to give the squared bicoherence
% bic2(f1,f2) = |B(f1,f2)|^2 / (<|X(f1)X(f2)|^2> <|X(f1+f2)|^2>), which is
% bounded in [0,1] by the Cauchy-Schwarz inequality.
%
% This is a different quantity from other 3rd-moment
% asymmetry/nonlinearity statistics (CO_trev, CO_tc3, SY_RampingWindows'
% asymAC1) which collapse all frequency structure into a single
% time-domain number per lag, so nonlinear coupling that's localized to a
% specific pair of frequency bands can average out to near zero. The
% bicoherence instead resolves quadratic phase coupling per frequency
% pair, directly detecting whether energy at f1 and f2 is
% phase-coupled to energy at f1+f2 (the frequency-domain signature of a
% quadratic nonlinearity). A single-lag time-domain moment
% cannot distinguish this from a nonlinearity spread evenly across all
% frequencies.
%
% ---INPUTS:
% y, the input time series
%
% segLength, the length (in samples) of each FFT segment. Segments
%            overlap by 50% and as many as fit in the series are
%            averaged over; the number of segments used, K, sets both
%            the variance of the bicoherence estimate and the
%            asymptotic significance threshold's precision (default: 64).
%
% maxN, the maximum number of samples to consider. Longer series are
%       cropped to their first maxN points, bounding runtime for
%       pathologically long series (cf. NL_RQA/NL_ZeroOneTest's maxN,
%       same rationale, though here cost scales ~linearly with N via the
%       number of segments K, rather than superlinearly, so this is a
%       safety margin rather than a necessity). Can be 'full' to disable
%       cropping (default).
%
% numSurr, the number of random-phase surrogates (SD_MakeSurrogates 'RP',
%          which preserve the power spectrum but destroy any phase
%          coupling) used to empirically calibrate the significance
%          threshold below, in place of its asymptotic approximation.
%          The threshold is the (1 - alpha) quantile of squared
%          bicoherence values pooled across all frequency pairs *and* all
%          surrogates -- pooling across the whole principal domain gives
%          a large effective sample even for a modest numSurr, since the
%          quantity being calibrated is a single global threshold, not a
%          per-triad one (default: 25).
%
% ---OUTPUTS:
% Mean, maximum, standard deviation, and skewness of the squared
% bicoherence over the non-redundant principal domain of frequency pairs
% (0 < f1 <= f2, f1 + f2 <= Nyquist); the mean squared bicoherence
% restricted to the self-coupling diagonal f1 = f2 (quadratic harmonic
% distortion -- the frequency-resolved analogue of CO_trev/CO_tc3's
% single time-domain statistics); the proportion of frequency pairs
% exceeding a surrogate-calibrated 95% significance threshold for
% quadratic phase coupling; the ratio of that empirical threshold to the
% standard analytic large-K approximation (K * bic2 ~ Exp(1)
% asymptotically under the null of a linear, ~Gaussian process, giving
% threshold -log(alpha)/K) -- a ratio far from 1 flags that the
% asymptotic approximation is untrustworthy for this series (e.g. because
% of non-stationarity), which real data was empirically found to trigger
% (see SP_Bicoherence's curation notes); and the Shannon entropy of the
% bicoherence surface (normalized to [0,1] by its maximum, i.e. the
% uniform-distribution entropy), reflecting whether any detected coupling
% is concentrated in a few frequency pairs or diffuse across many (cf.
% "bispectral index"-type measures used in EEG analysis).

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
%% Check inputs, set defaults:
% ------------------------------------------------------------------------------
if size(y, 2) > size(y, 1)
    y = y'; % Time series must be a column vector
end

if nargin < 2 || isempty(segLength)
    segLength = 64;
end
if nargin < 3 || isempty(maxN)
    maxN = 'full';
end
if nargin < 4 || isempty(numSurr)
    numSurr = 25;
end

minSegLength = 16;
if segLength < minSegLength
    error('segLength = %u is too short for a meaningful FFT segment (need >= %u)', segLength, minSegLength);
end

N = length(y);
if ischar(maxN) && strcmp(maxN, 'full')
    % No cropping
elseif N > maxN
    warning('Time series (%u samples) exceeds maxN = %u; analyzing the first %u samples', N, maxN, maxN);
    y = y(1:maxN);
    N = maxN;
end

% ------------------------------------------------------------------------------
%% Segment geometry (50% overlap, as many segments as fit) -- shared by the
%% real series and every surrogate, since they're all the same length
% ------------------------------------------------------------------------------
step = floor(segLength / 2);
minNumSeg = 8; % need enough segments for a meaningful bicoherence estimate
numSeg = floor((N - segLength) / step) + 1;

if numSeg < minNumSeg
    warning(['Time series (N = %u) too short for segLength = %u to form >= %u ' ...
             '50%%-overlapping segments'], N, segLength, minNumSeg);
    out = NaN; return
end

halfN = floor(segLength / 2) + 1; % index i (1..halfN) <-> frequency (i-1)/segLength, up to Nyquist
win = hamming(segLength);

[I, J] = ndgrid(1:halfN, 1:halfN); % I(i,j) = i, J(i,j) = j
sumIdx = I + J - 1; % index of f1 + f2 (bin arithmetic: index (i-1)+(j-1)+1)
validSum = (sumIdx <= halfN); % triads whose sum frequency is within Nyquist
mask = validSum & (J >= I) & (I >= 2); % non-redundant principal domain, excluding DC

geom = struct('step', step, 'numSeg', numSeg, 'segLength', segLength, ...
              'halfN', halfN, 'win', win, 'I', I, 'J', J, 'sumIdx', sumIdx, 'validSum', validSum);

% ------------------------------------------------------------------------------
%% Bicoherence of the real series
% ------------------------------------------------------------------------------
bic2 = SUB_BicoherenceGrid(y, geom);
bicVals = bic2(mask);

if isempty(bicVals) || all(~isfinite(bicVals))
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Summary statistics
% ------------------------------------------------------------------------------
out.meanBic = mean(bicVals);
out.maxBic = max(bicVals);
out.stdBic = std(bicVals);
out.skewBic = skewness(bicVals);

% Normalized Shannon entropy of the bicoherence surface (0 = all coupling
% concentrated in one frequency pair, 1 = uniformly diffuse):
p = bicVals(bicVals > 0);
p = p / sum(p);
out.entropy = -sum(p .* log(p)) / log(length(bicVals));

% Self-coupling diagonal (f1 = f2, quadratic harmonic distortion):
diagMask = validSum & (I == J) & (I >= 2);
diagVals = bic2(diagMask);
out.meanBicDiag = mean(diagVals);

% ------------------------------------------------------------------------------
%% Surrogate-calibrated significance threshold
% ------------------------------------------------------------------------------
% Random-phase surrogates preserve the power spectrum (linear structure)
% but destroy quadratic phase coupling -- exactly the null hypothesis a
% bicoherence significance test needs. Pool their bic2 values across the
% whole principal domain to estimate the alpha = 0.05 threshold
% empirically, rather than relying on the analytic large-K approximation
% alone (which assumes segment-to-segment i.i.d. complex Gaussian
% behavior under the null -- an assumption real, non-stationary data can
% violate).
alpha = 0.05;
surrogates = SD_MakeSurrogates(y, 'RP', numSurr, [], 'default');

nullVals = cell(numSurr, 1);
for s = 1:numSurr
    bic2surr = SUB_BicoherenceGrid(surrogates(:, s), geom);
    nullVals{s} = bic2surr(mask);
end
nullVals = vertcat(nullVals{:});
nullVals = nullVals(isfinite(nullVals));

surrThresh = quantile(nullVals, 1 - alpha);
out.propSig = mean(bicVals > surrThresh);

% How far the standard asymptotic threshold (-log(alpha)/K, from
% K * bic2 ~ Exp(1) under the null) is from the empirical one -- a ratio
% far from 1 means the asymptotic approximation shouldn't be trusted for
% this series:
analyticThresh = -log(alpha) / numSeg;
out.threshRatio = surrThresh / analyticThresh;

end

% ------------------------------------------------------------------------------
function bic2 = SUB_BicoherenceGrid(y, geom)
    % Segment-averaged squared bicoherence grid for series y, given the
    % (series-length-independent-in-shape) segmentation geometry geom.
    halfN = geom.halfN;
    Bnum = zeros(halfN, halfN); % complex triple-product sum
    P12 = zeros(halfN, halfN); % sum |X(f1) X(f2)|^2
    P3 = zeros(halfN, halfN); % sum |X(f1+f2)|^2

    for k = 1:geom.numSeg
        startIdx = (k - 1) * geom.step + 1;
        seg = y(startIdx:startIdx + geom.segLength - 1);
        seg = seg - mean(seg); % demean each segment before windowing
        X = fft(seg .* geom.win, geom.segLength);
        Xh = X(1:halfN); % one-sided spectrum, DC to Nyquist

        outerXX = Xh(geom.I) .* Xh(geom.J); % X(f1) * X(f2) for every (f1,f2) pair
        Xsum = zeros(halfN, halfN);
        Xsum(geom.validSum) = Xh(geom.sumIdx(geom.validSum)); % X(f1+f2) where defined

        triple = outerXX .* conj(Xsum);
        triple(~geom.validSum) = 0;

        Bnum = Bnum + triple;
        P12 = P12 + abs(outerXX).^2;
        P3(geom.validSum) = P3(geom.validSum) + abs(Xh(geom.sumIdx(geom.validSum))).^2;
    end

    bic2 = abs(Bnum).^2 ./ (P12 .* P3 + eps); % squared bicoherence, bounded in [0,1]
end
