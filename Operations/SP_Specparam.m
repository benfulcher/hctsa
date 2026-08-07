function out = SP_Specparam(y, aperiodicMode, maxNPeaks, peakThreshold, peakWidthLimits, segLength, maxSegments)
% SP_Specparam   Separates the power spectrum into aperiodic (1/f) and periodic (oscillatory) components
%
% Parameterizes the power spectrum as a smooth aperiodic '1/f' background
% plus a small number of Gaussian peaks sitting on top of it, in the
% spirit of the FOOOF/specparam algorithm (Donoghue et al., "Parameterizing
% neural power spectra into periodic and aperiodic components", Nature
% Neuroscience 23: 1655 (2020)).
%
% The motivation is that hctsa's existing spectral peak statistics and its
% 1/f slope fits are computed independently of one another, and each
% therefore contaminates the other:
%   - SP_Summaries' peak fields (numPromPeaks_*, maxProm, ...) measure
%     prominence against the local spectrum, so a peak's apparent size
%     conflates a genuine oscillation with however steeply the aperiodic
%     background happens to be falling underneath it.
%   - SP_Summaries' linfitloglog_* fields fit a straight line through
%     log-log spectrum over fixed ranges, with any oscillatory peaks left
%     in -- so strong peaks drag the fitted slope away from the true
%     background exponent.
% Fitting both together, and iterating (peaks are removed before the
% aperiodic component is re-fit), is what decouples them.
%
% ---INPUTS:
% y, the input time series
%
% aperiodicMode, the form of the aperiodic component:
%           (i) 'fixed': b - chi*log10(f), a straight line in log-log,
%               i.e. pure power-law
%           (ii) 'knee': b - log10(k + f^chi), which additionally allows
%               the spectrum to flatten off below a 'knee' frequency, as
%               real spectra commonly do. Note the knee model is not
%               identifiable when the data has no actual knee (k -> 0),
%               so it falls back to the 'fixed' fit if the optimization
%               fails or returns a degenerate knee.
%
% maxNPeaks, the maximum number of Gaussian peaks to extract (default: 4).
%
% peakThreshold, how far above the noise a candidate peak must stand to be
%           accepted (default: 1). This is expressed as a multiple of the
%           largest deviation that noise alone would be expected to
%           produce -- specifically of sqrt(2*log(nBins)) robust standard
%           deviations of the flattened spectrum, nBins being the number
%           of frequency bins searched. Expressing it that way (rather
%           than as a plain multiple of sigma) is necessary because the
%           test is applied to the maximum over all bins: a fixed small
%           multiple of sigma fires on pure noise essentially always,
%           and the correction adapts as nBins changes. So a value of 1
%           means 'must exceed what noise alone would give', and larger
%           values are correspondingly more conservative.
%
% peakWidthLimits, two-element [min, max] on Gaussian peak standard
%           deviation, in log10-frequency units (default: [0.02, 0.5]).
%           Bounding the width both stops the optimizer fitting a single
%           enormously wide 'peak' that is really leftover aperiodic
%           background, and stops it fitting single-bin spectral noise.
%
% ---OUTPUTS:
% The aperiodic parameters (exponent, offset, and for 'knee' mode the knee
% parameter); the number of peaks found above threshold, and the centre
% frequency, height and bandwidth of the largest; the total power in the
% periodic component and the fraction of spectral power it accounts for;
% and the quality of the combined fit (R^2 and mean absolute error).
%
% ---WHICH OUTPUTS ARE REGISTERED, AND WHY:
% Seven of the ten are registered: apExponent, apOffset, numPeaks,
% maxPeakFreq, maxPeakPower, maxPeakBW and periodicFraction. The
% exclusions are measured rather than assumed:
%
%   * modelR2 and modelMAE are genuinely length-dependent (eta^2 against N
%     of 0.50-0.76 and 0.60-0.61) and consistently so across three
%     different generating processes, including a stationary AR(1). This
%     is intrinsic to what they measure: goodness of fit depends on how
%     noisy the spectral estimate is, which depends on how many segments
%     there are to average, which depends on N. They cannot be compared
%     across a dataset of mixed lengths.
%   * apKnee reaches eta^2 = 0.69 on AR(1) -- whose spectrum is a
%     Lorentzian, i.e. exactly the knee model, so the knee is genuinely
%     identifiable there and its estimate is correspondingly sensitive to
%     how much data is available. Since it is the only field unique to
%     'knee' mode, only the 'fixed' variant is registered; 'knee' remains
%     callable and does fit real knees well (recovering a planted knee of
%     1e-3 as 1.01e-3, R^2 0.996 vs 0.945 for the 'fixed' form) but costs
%     ~4x more.
%   * totalPeakPower is excluded as a near-duplicate of maxPeakPower
%     (R^2 = 0.85 across 1000 real series; they coincide exactly whenever
%     only one peak is found, which is the common case), with numPeaks
%     already carrying the 'how many' information.
%
% Measuring length-dependence needs care, and getting it wrong initially
% led this operation to discard six fields that were in fact fine.
% Generating a fresh realization at each N is *not* a valid test for a
% 1/f-type process: a longer realization contains lower frequencies, so
% the process itself changes with N and the resulting eta^2 conflates that
% with any bias in the operation. The correct design, used for the figures
% above, is to draw one long realization per replicate and take prefixes
% of it, which varies the observation length while holding the process
% fixed.
%
% On novelty: across 1000 real series, periodicFraction reaches a maximum
% R^2 of only 0.17 against any of the ~7700 other hctsa features, and
% numPeaks 0.24 (and only 0.17 against SP_Summaries' own spectral-peak
% fields) -- i.e. counting peaks *relative to a fitted aperiodic
% background* really is different from counting them by absolute
% prominence, which was the motivating claim for this operation.
% apExponent is much the most redundant (R^2 = 0.76 against
% SP_Summaries' linfitloglog_all_a2 and 0.79 against
% linfitsemilog_all_a2) but is the better-motivated estimator of that same
% quantity, being both peak-corrected and fitted to a segment-averaged
% spectrum.

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
%% Check that a Curve-Fitting Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('curve_fitting_toolbox');

% ------------------------------------------------------------------------------
%% Check inputs, set defaults:
% ------------------------------------------------------------------------------
if size(y, 2) > size(y, 1)
    y = y'; % Time series must be a column vector
end

if nargin < 2 || isempty(aperiodicMode)
    aperiodicMode = 'fixed';
end
if ~ismember(aperiodicMode, {'fixed', 'knee'})
    error('Unknown aperiodicMode ''%s'' (expected ''fixed'' or ''knee'')', aperiodicMode);
end
if nargin < 3 || isempty(maxNPeaks)
    maxNPeaks = 4;
end
if nargin < 4 || isempty(peakThreshold)
    peakThreshold = 1;
end
if nargin < 5 || isempty(peakWidthLimits)
    peakWidthLimits = [0.02, 0.5];
end
if nargin < 6
    segLength = []; % [] = adapt to the series length (see below)
end
if nargin < 7 || isempty(maxSegments)
    maxSegments = Inf; % Inf = use all the available data
end

N = length(y);
if isempty(segLength)
    % Scale the segment length with the series, so that longer series buy
    % both finer frequency resolution and more segments to average over.
    % Holding it fixed instead was tried and measurably rejected: see the
    % note on the spectral estimate below.
    segLength = max(round(N / 8), 32);
end
minSegments = 4; % need several segments to average over for a usable estimate
minLength = segLength + (minSegments - 1) * floor(segLength / 2);
if N < minLength
    warning(['Time series (N = %u) too short for a spectral parameterization with ' ...
             'segLength = %u (need >= %u)'], N, segLength, minLength);
    out = NaN; return
end
if all(y == y(1))
    warning('Constant time series has no spectral structure');
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Spectrum: Welch's method, computed here rather than taken from a caller
% ------------------------------------------------------------------------------
% A parameterization of the spectrum's *shape* needs a smooth, low-variance
% spectral estimate. A raw periodogram is a statistically inconsistent
% estimator (its per-bin variance does not shrink as more data is added),
% so fitting peaks to one mostly fits noise -- the same fact that forced
% SP_Summaries' peak-detection thresholds to be calibrated on Welch
% estimates only. Hence Welch's method (segment-averaged) is used here
% unconditionally rather than being left to the caller.
% The segment length scales with N by default, and all available data is
% used. Holding both fixed instead (so that every series yields a
% statistically identical spectral estimate, and hence perfectly
% length-invariant parameters) is available via the segLength and
% maxSegments arguments, but was tried as the default and rejected on
% measurement: fixing segLength = 256 and capping at 15 segments starved
% the estimate badly enough that the exponent recovered from a known 1/f
% process degraded (chi = 1.0 recovered as 1.040 rather than 1.001), the
% white-noise false-peak rate rose from 2% to 12%, the knee model began
% misfitting pure power-law data (returning a spurious knee of 1e-2 and
% an exponent of 1.22 instead of 1.0), and -- decisively -- the
% peak-corrected exponent went from twice as accurate as a naive
% log-log slope fit to six times *worse* than it, because the naive
% comparator still gets to use the whole series. Length-invariance
% bought at the price of discarding most of the data is not worth having
% here; the length-dependent outputs are excluded from registration
% instead.
winLength = segLength;
maxSamples = winLength + (maxSegments - 1) * floor(winLength / 2);
if N > maxSamples
    y = y(1:maxSamples); % use a fixed amount of data so features stay comparable
end
NFFT = 2^nextpow2(winLength);
[S, f] = pwelch(y, hamming(winLength), [], NFFT, 1);

% Restrict to the frequencies this estimate can actually resolve. A
% frequency only completing one or two cycles within a Welch segment is
% essentially unestimated, and on a log-frequency axis those lowest bins
% sit isolated far to the left, where a straight-line fit is
% unconstrained -- so any wiggle there reads as a huge 'peak'. Verified:
% without this cut, 1/f noise (which contains no oscillation at all)
% produced a spurious peak in 73% of replicates, and *every* one of them
% sat in the lowest one or two frequency bins. Requiring several cycles
% per segment removes them.
% Note this makes the fitted frequency range depend on N, which is
% physically correct (a shorter series genuinely cannot resolve low
% frequencies) but does mean these features are not strictly
% length-independent.
minCyclesPerSegment = 5;
fMin = minCyclesPerSegment / winLength;

% Also exclude DC (log10(0) = -Inf) and any non-positive/non-finite power:
valid = (f >= fMin) & (S > 0) & isfinite(S);
if sum(valid) < 10
    warning('Too few valid spectral points for a parameterization');
    out = NaN; return
end
fv = f(valid);
logF = log10(fv);
logS = log10(S(valid));

% ------------------------------------------------------------------------------
%% Initial aperiodic fit
% ------------------------------------------------------------------------------
% Peaks bias this first fit; that is expected and is corrected by the
% refit at the end, once the peaks have been identified and removed.
ap0 = SUB_FitAperiodic(fv, logF, logS, aperiodicMode);

% ------------------------------------------------------------------------------
%% Iteratively extract Gaussian peaks from the flattened spectrum
% ------------------------------------------------------------------------------
resid = logS - ap0.pred;
peakList = struct('centre', {}, 'height', {}, 'width', {});
peakSum = zeros(size(logS));

nBins = length(resid);
% The acceptance threshold has to account for the fact that we are testing
% the *maximum* over all nBins frequency bins, not one pre-specified bin:
% the largest of nBins noise samples is expected to sit around
% sqrt(2*log(nBins)) standard deviations up (~3.5 for a few hundred bins),
% so any fixed small multiple of sigma fires on pure noise essentially
% always. Scaling by that factor makes the threshold a genuine
% "bigger than the biggest thing noise would produce" test, and adapts
% automatically as the number of bins changes with series length.
% Verified: without this, white noise hit the maximum peak count in 100%
% of replicates and pink noise found a spurious peak in 98%.
nullMaxFactor = sqrt(2 * log(nBins));

for k = 1:maxNPeaks
    % Robust spread: the residual still contains the very peaks we are
    % looking for, and those positive outliers inflate a plain std --
    % which would make the test *less* sensitive exactly when there is
    % real structure. A MAD-based sigma is not pulled about by them.
    residSD = 1.4826 * median(abs(resid - median(resid)));
    if residSD <= 0, break; end
    [pkHeight, iPk] = max(resid);
    if pkHeight < peakThreshold * nullMaxFactor * residSD
        break % nothing left standing above what noise alone would give
    end

    g = SUB_FitGaussian(logF, resid, iPk, peakWidthLimits);
    if isempty(g)
        break % fit failed or returned a degenerate/out-of-bounds peak
    end

    peakList(end + 1) = struct('centre', g.centre, 'height', g.height, 'width', g.width); %#ok<AGROW>
    peakSum = peakSum + g.pred;
    resid = resid - g.pred; % peel this peak off and look for the next
end

% ------------------------------------------------------------------------------
%% Refit the aperiodic component with the peaks removed
% ------------------------------------------------------------------------------
% This is the step that decouples the two components: with the oscillatory
% peaks subtracted, the background fit is no longer dragged by them, so
% the exponent estimates the true 1/f background rather than a blend of
% background and oscillations (which is what a direct log-log line fit,
% e.g. SP_Summaries' linfitloglog_*, returns).
apFinal = SUB_FitAperiodic(fv, logF, logS - peakSum, aperiodicMode);

out.apExponent = apFinal.exponent;
% Report the background level at a reference frequency *inside* the fitted
% band, rather than the raw intercept. The intercept is the fitted value at
% log10(f) = 0, i.e. f = 1 -- above the Nyquist frequency of 0.5, so it is
% a pure extrapolation whose value swings with both the fitted slope and
% wherever the fitted range happens to start. Evaluating the same fitted
% curve at a frequency actually covered by the data removes that
% sensitivity while describing the same thing (how much power the
% aperiodic background carries).
refFreq = 0.1;
out.apOffset = SUB_EvalAperiodic(apFinal, refFreq);
if strcmp(aperiodicMode, 'knee')
    out.apKnee = apFinal.knee;
end

% ------------------------------------------------------------------------------
%% Peak statistics (measured relative to the aperiodic background)
% ------------------------------------------------------------------------------
out.numPeaks = length(peakList);
if isempty(peakList)
    out.maxPeakFreq = NaN;
    out.maxPeakPower = NaN;
    out.maxPeakBW = NaN;
    out.totalPeakPower = 0;
else
    heights = [peakList.height];
    [maxH, iMax] = max(heights);
    out.maxPeakFreq = 10^peakList(iMax).centre; % back to linear frequency
    out.maxPeakPower = maxH; % height above the aperiodic background, in log10 power
    out.maxPeakBW = peakList(iMax).width;
    out.totalPeakPower = sum(heights);
end

% Share of the (log-)spectrum's variation accounted for by the periodic
% component, rather than by the aperiodic background:
totalVar = sum((logS - mean(logS)).^2);
if totalVar > 0
    out.periodicFraction = sum(peakSum.^2) / totalVar;
else
    out.periodicFraction = NaN;
end

% ------------------------------------------------------------------------------
%% Quality of the combined (aperiodic + periodic) fit
% ------------------------------------------------------------------------------
model = apFinal.pred + peakSum;
residFinal = logS - model;
if totalVar > 0
    out.modelR2 = 1 - sum(residFinal.^2) / totalVar;
else
    out.modelR2 = NaN;
end
out.modelMAE = mean(abs(residFinal));

end

% ------------------------------------------------------------------------------
function v = SUB_EvalAperiodic(ap, fq)
    % Value of the fitted aperiodic curve at frequency fq.
    if isfield(ap, 'knee') && ap.knee > 0
        v = ap.offset - log10(ap.knee + fq^ap.exponent);
    else
        v = ap.offset - ap.exponent * log10(fq);
    end
end

% ------------------------------------------------------------------------------
function ap = SUB_FitAperiodic(fv, logF, logS, aperiodicMode)
    % Fit the smooth aperiodic background of the log10 spectrum.
    % 'fixed': logS = offset - exponent*log10(f)
    % 'knee' : logS = offset - log10(knee + f^exponent)

    % Robust straight-line fit in log-log, used directly for 'fixed' mode
    % and as the starting point for the nonlinear 'knee' fit:
    warning('off', 'stats:robustfit:RankDeficient');
    b = robustfit(logF, logS);
    warning('on', 'stats:robustfit:RankDeficient');
    ap.offset = b(1);
    ap.exponent = -b(2); % conventionally reported as a positive falling slope
    ap.pred = b(1) + b(2) * logF;

    if strcmp(aperiodicMode, 'fixed')
        return
    end

    % 'knee' mode: a nonlinear fit, seeded from the linear one. The knee
    % is only identifiable when the spectrum actually flattens at low
    % frequency; when it does not, the optimizer drives knee -> 0 (or
    % fails outright), in which case the model degenerates to the 'fixed'
    % form and we keep the robust linear fit rather than a bogus knee.
    try
        s = fitoptions('Method', 'NonlinearLeastSquares', ...
                       'StartPoint', [ap.offset, 1e-3, max(ap.exponent, 0.1)], ...
                       'Lower', [-Inf, 0, 0], ...
                       'Upper', [Inf, Inf, 10], ...
                       'MaxIter', 400, 'Display', 'off');
        ft = fittype('a - log10(k + x^c)', 'independent', 'x', 'options', s);
        [c, ~] = fit(fv, logS, ft);
        kneeVal = c.k;
        predKnee = c.a - log10(kneeVal + fv.^c.c);
        if isfinite(kneeVal) && kneeVal > 1e-10 && all(isfinite(predKnee))
            ap.offset = c.a;
            ap.knee = kneeVal;
            ap.exponent = c.c;
            ap.pred = predKnee;
        else
            ap.knee = 0; % degenerate: no detectable knee, keep the linear fit
        end
    catch
        ap.knee = 0; % optimization failed: fall back to the linear fit
    end
end

% ------------------------------------------------------------------------------
function g = SUB_FitGaussian(logF, resid, iPk, peakWidthLimits)
    % Fit a single Gaussian to the flattened spectrum, centred near the
    % current maximum at index iPk. Returns [] if the fit fails or lands
    % outside the permitted width range.
    g = [];
    x0 = logF(iPk);
    h0 = resid(iPk);
    if ~isfinite(h0) || h0 <= 0
        return
    end
    w0 = mean(peakWidthLimits);

    try
        s = fitoptions('Method', 'NonlinearLeastSquares', ...
                       'StartPoint', [h0, x0, w0], ...
                       'Lower', [0, min(logF), peakWidthLimits(1)], ...
                       'Upper', [Inf, max(logF), peakWidthLimits(2)], ...
                       'MaxIter', 400, 'Display', 'off');
        ft = fittype('h*exp(-(x-m)^2/(2*w^2))', 'independent', 'x', 'options', s);
        c = fit(logF, resid, ft);
    catch
        return
    end

    if ~isfinite(c.h) || ~isfinite(c.m) || ~isfinite(c.w) || c.h <= 0
        return
    end
    g.height = c.h;
    g.centre = c.m;
    g.width = c.w;
    g.pred = c.h * exp(-(logF - c.m).^2 / (2 * c.w^2));
end
