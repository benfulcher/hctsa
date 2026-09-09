function out = SP_Rhythmicity(y, m, tau, nFreq, fLow, fHigh, smoothWin)
% SP_Rhythmicity   Phase-persistence ("rhythmicity") spectrum of a time series.
%
% A power spectrum cannot tell a sustained oscillation apart from a train of
% brief, high-power bursts: both deposit energy at the same frequency. This
% operation instead asks, at each frequency, how *predictable* the phase is a
% fixed number of cycles into the future -- high for an oscillation that keeps
% ticking (sustained/"rhythmic"), low for one that only flickers on
% intermittently (transient/"bursty") even when the two have identical power.
%
% Method (adapted from G. Karvat et al., "Universal rhythmic architecture
% uncovers two modes of neural dynamics", Nature Communications (2026), the
% Lagged-Angle Vector Index / LAVI; itself a single-lag refinement of the
% lagged-coherence measure of J. Fransen et al., NeuroImage 111, 158 (2015)):
%
% (1) Complex Morlet wavelet transform at each analysis frequency f (m cycles
%     per wavelet), giving instantaneous phase phi(t,f) and amplitude a(t,f).
%     No Wavelet Toolbox dependency -- the transform is an FFT-domain
%     convolution here.
% (2) At each f, compare the wavelet coefficient at time t with the one at
%     t+L, where the lag L = round(tau/f) samples is a fixed tau cycles. Two
%     rhythmicity spectra are built from this comparison:
%       * LAVI(f) = | mean_t exp( 1i*[ phi(t,f) - phi(t+L,f) ] ) |, the mean
%         resultant length of the phase-difference distribution (pure phase,
%         in [0,1]); 1 => phase advances by a perfectly repeatable amount
%         every tau cycles (a stable rhythm), ~0 => the phase L samples ahead
%         is unpredictable.
%       * LC(f), the amplitude-weighted version: the classic lagged-coherence
%         estimator | sum_t x(t,f) conj(x(t+L,f)) | / sqrt( sum|x(t,f)|^2 *
%         sum|x(t+L,f)|^2 ), also in [0,1]. LAVI weights every timepoint
%         equally; LC lets high-amplitude timepoints dominate, so the two
%         diverge when rhythmic and arrhythmic epochs differ in power.
% (3) Each spectrum has a stable baseline (its own median across frequency),
%     so peaks above baseline flag sustained bands and troughs below it flag
%     transient/bursty bands. Both are summarized below with the same battery
%     of statistics (prefixes lavi_ and lc_), so the phase-only and
%     amplitude-weighted views can be compared directly.
%
% Following hctsa convention the series is treated as evenly sampled with unit
% spacing: frequencies are in cycles per sample (0 to 0.5) and reported
% angular frequencies are 2*pi*f, as in SP_Summaries. A given frequency is
% skipped (NaN, excluded from the summaries) when the series is too short to
% support its wavelet plus lag; the operation returns a single NaN if fewer
% than 10 frequencies survive.
%
% cf. SP_SpectralTimeFreq, whose spectral kurtosis also targets bursty
% frequency bands but via the variance of band power across a spectrogram's
% windows -- a power-domain construct that says nothing about phase
% continuity; and SP_SummariesPhase, which characterizes the static Fourier
% phase spectrum of the whole series rather than phase persistence over time.
%
% ---INPUTS:
% y, the input time series
%
% m, width of the Morlet wavelet in cycles (default: 5). Larger m = sharper
%    frequency resolution, blunter time resolution.
%
% tau, the lag at which phase persistence is measured, in cycles (default:
%      1.5). Per Karvat et al. the relative shape of the rhythmicity spectrum
%      is nearly independent of this choice, so a single lag suffices.
%
% nFreq, number of log-spaced analysis frequencies (default: 40).
%
% fLow, lowest analysis frequency in cycles/sample (default: 0.02).
%
% fHigh, highest analysis frequency in cycles/sample (default: 0.4; must be
%        below the Nyquist frequency of 0.5).
%
% smoothWin, width (in frequency bins, forced odd) of the moving-average
%       (Daniell) spectral window applied before counting rhythmic bands
%       (default: 5; 1 disables it). Welch-style segment averaging is not
%       used -- for typical hctsa series lengths each segment would fall
%       below the wavelet-plus-lag length -- so estimate variance is instead
%       reduced by this frequency-domain smoother, and only the band count
%       (lavi_numPeaks / lc_numPeaks) is taken from the smoothed spectrum;
%       every other statistic uses the raw spectrum.
%
% ---OUTPUTS: the same statistics are returned for the phase-only spectrum
% LAVI(f) (prefix lavi_) and the amplitude-weighted spectrum LC(f) (prefix
% lc_), computed over the frequencies that survive the length check:
% *_mean, *_median, *_std, *_iqr, *_max, *_min, *_range, basic level and
%       spread of rhythmicity across frequency (*_median is the baseline).
% *_maxOnMedian, peak rhythmicity relative to the baseline -- how far the
%       most sustained band stands out.
% *_freqAtMax, *_freqAtMin, angular frequency (2*pi*f) of the most sustained
%       and most transient band.
% *_meanAboveMedian, mean amount by which above-baseline frequencies exceed
%       the baseline (average prominence of sustained bands).
% *_fracSustained, *_fracTransient, fraction of frequencies sitting more than
%       2 scaled-MAD above / below the baseline -- a data-driven stand-in for
%       the surrogate "noise ribbon" of the original method (which needs many
%       minutes of data and is not attempted here).
% *_areaSustained, area of the spectrum above its baseline, integrated over
%       log-frequency and normalized by the log-frequency span: the total
%       "sustained mass" of the spectrum.
% *_slope, ordinary-least-squares slope of the spectrum against log10(f):
%       does rhythmicity concentrate at low or high frequencies?
% *_ac1, lag-1 autocorrelation of the spectrum across frequency: a smooth
%       profile (broad bands) vs. a jagged one.
% *_numPeaks, number of local maxima of the smoothed spectrum that exceed the
%       baseline by at least one scaled-MAD (number of distinct rhythmic
%       bands).
% *_corrPower, Spearman correlation between the spectrum and log10 of the
%       mean wavelet power at f: near 0 => rhythmicity carries information
%       beyond the power spectrum; near 1 => it largely tracks power.
% lc_corrLAVI, Spearman correlation between the LC and LAVI spectra: near 1 =>
%       amplitude weighting barely changes the picture; lower => rhythmic and
%       arrhythmic epochs differ systematically in power.

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
if nargin < 2 || isempty(m)
    m = 5;
end
if nargin < 3 || isempty(tau)
    tau = 1.5;
end
if nargin < 4 || isempty(nFreq)
    nFreq = 40;
end
if nargin < 5 || isempty(fLow)
    fLow = 0.02;
end
if nargin < 6 || isempty(fHigh)
    fHigh = 0.4;
end
if nargin < 7 || isempty(smoothWin)
    smoothWin = 5;
end
smoothWin = max(1, round(smoothWin));
if mod(smoothWin, 2) == 0
    smoothWin = smoothWin + 1; % keep the moving average centered
end

if fHigh >= 0.5
    error('fHigh = %.3f must be below the Nyquist frequency (0.5)', fHigh);
end
if fLow <= 0 || fLow >= fHigh
    error('Need 0 < fLow < fHigh');
end

N = length(y);
if std(y) < eps
    out = NaN; return % constant series carries no phase information
end
y = y - mean(y);

freqs = logspace(log10(fLow), log10(fHigh), nFreq)';

% ------------------------------------------------------------------------------
%% Complex Morlet wavelet transform (FFT-domain convolution; no toolbox)
% ------------------------------------------------------------------------------
% Per frequency we also record the wavelet half-length, so the lagged step
% below can discard the edge-contaminated samples.
X = complex(nan(N, nFreq));
halfLen = nan(nFreq, 1);
for fi = 1:nFreq
    f = freqs(fi);
    sigmaT = m / (2*pi*f); % Gaussian envelope s.d. in samples
    h = ceil(3.5 * sigmaT); % wavelet support: +/- 3.5 s.d.
    if 2*h + 1 >= N
        continue % wavelet longer than the series: leave this frequency NaN
    end
    halfLen(fi) = h;

    t = (-h:h)';
    psi = exp(2i*pi*f*t) .* exp(-t.^2 / (2*sigmaT^2));
    psi = psi - mean(psi); % zero-mean (wavelet admissibility)
    psi = psi / sqrt(sum(abs(psi).^2)); % unit energy

    % Linear convolution via FFT, then take the length-N 'same' segment:
    Lconv = N + numel(psi) - 1;
    c = ifft(fft(y, Lconv) .* fft(psi, Lconv));
    X(:, fi) = c(h+1 : h+N);
end

% ------------------------------------------------------------------------------
%% Lagged phase persistence at each frequency: LAVI (phase-only) and LC
%% (amplitude-weighted lagged coherence)
% ------------------------------------------------------------------------------
lavi = nan(nFreq, 1);
lc = nan(nFreq, 1);
logPower = nan(nFreq, 1);
for fi = 1:nFreq
    if isnan(halfLen(fi))
        continue
    end
    f = freqs(fi);
    h = halfLen(fi);
    L = max(1, round(tau / f)); % lag in samples

    lo = h + 1;
    hi = N - h - L;
    minValid = max(64, 3*L); % enough independent lag spans for a stable estimate
    if hi - lo + 1 < minValid
        continue
    end

    x1 = X(lo:hi, fi);
    x2 = X(lo+L:hi+L, fi);

    lavi(fi) = abs(mean(exp(1i * (angle(x1) - angle(x2)))));
    lc(fi) = abs(sum(x1 .* conj(x2))) / sqrt(sum(abs(x1).^2) * sum(abs(x2).^2));
    logPower(fi) = log10(mean(abs(x1).^2));
end

% ------------------------------------------------------------------------------
%% Summarize each rhythmicity spectrum
% ------------------------------------------------------------------------------
valid = isfinite(lavi) & isfinite(lc);
if sum(valid) < 10
    out = NaN; return % too little of the spectrum is estimable to summarize
end

f = freqs(valid);
logf = log10(f);
P = logPower(valid);
Lavi = lavi(valid);
Lc = lc(valid);

out = SUB_summarize(Lavi, f, logf, P, smoothWin, 'lavi_');
out = SUB_mergeStruct(out, SUB_summarize(Lc, f, logf, P, smoothWin, 'lc_'));

% How much does the amplitude weighting change the rhythmicity spectrum?
if std(Lavi) < eps || std(Lc) < eps
    out.lc_corrLAVI = NaN;
else
    out.lc_corrLAVI = SUB_spearman(Lc, Lavi);
end

end

% ==============================================================================
function s = SUB_summarize(S, f, logf, P, smoothWin, prefix)
    % Battery of summary statistics for one rhythmicity spectrum S(f).
    med = median(S);
    scaledMAD = 1.4826 * median(abs(S - med));

    s.([prefix 'mean']) = mean(S);
    s.([prefix 'median']) = med;
    s.([prefix 'std']) = std(S);
    s.([prefix 'iqr']) = iqr(S);
    s.([prefix 'max']) = max(S);
    s.([prefix 'min']) = min(S);
    s.([prefix 'range']) = max(S) - min(S);
    s.([prefix 'maxOnMedian']) = max(S) / med;

    [~, iMax] = max(S);
    [~, iMin] = min(S);
    s.([prefix 'freqAtMax']) = 2*pi*f(iMax);
    s.([prefix 'freqAtMin']) = 2*pi*f(iMin);

    above = S(S > med) - med;
    if isempty(above)
        s.([prefix 'meanAboveMedian']) = 0;
    else
        s.([prefix 'meanAboveMedian']) = mean(above);
    end

    if scaledMAD < eps
        s.([prefix 'fracSustained']) = 0;
        s.([prefix 'fracTransient']) = 0;
    else
        s.([prefix 'fracSustained']) = mean(S > med + 2*scaledMAD);
        s.([prefix 'fracTransient']) = mean(S < med - 2*scaledMAD);
    end

    % Area above baseline, integrated over log-frequency, span-normalized:
    s.([prefix 'areaSustained']) = trapz(logf, max(S - med, 0)) / (logf(end) - logf(1));

    % OLS slope of the spectrum against log-frequency:
    pf = polyfit(logf, S, 1);
    s.([prefix 'slope']) = pf(1);

    % Lag-1 autocorrelation of the spectrum across frequency:
    S0 = S - mean(S);
    denom = sum(S0.^2);
    if denom < eps
        s.([prefix 'ac1']) = NaN;
    else
        s.([prefix 'ac1']) = sum(S0(1:end-1) .* S0(2:end)) / denom;
    end

    % Number of rhythmic bands: local maxima of the smoothed spectrum that
    % clear the baseline by at least one scaled-MAD (a Daniell/moving-average
    % window stands in for Welch-style segment averaging, which the series is
    % too short to afford here).
    Ssm = movmean(S, smoothWin);
    thresh = med + max(scaledMAD, 0);
    dS = diff(Ssm);
    isPeak = [false; (dS(1:end-1) > 0) & (dS(2:end) <= 0); false];
    s.([prefix 'numPeaks']) = sum(isPeak & (Ssm > thresh));

    % Does rhythmicity just track power?
    if std(P) < eps
        s.([prefix 'corrPower']) = NaN;
    else
        s.([prefix 'corrPower']) = SUB_spearman(S, P);
    end
end

% ==============================================================================
function a = SUB_mergeStruct(a, b)
    fn = fieldnames(b);
    for i = 1:numel(fn)
        a.(fn{i}) = b.(fn{i});
    end
end

% ==============================================================================
function r = SUB_spearman(a, b)
    % Spearman rank correlation without a Statistics Toolbox dependency
    % (continuous wavelet quantities: ties are measure-zero, ignored).
    n = numel(a);
    [~, oa] = sort(a); ra = zeros(n,1); ra(oa) = 1:n;
    [~, ob] = sort(b); rb = zeros(n,1); rb(ob) = 1:n;
    ra = ra - mean(ra);
    rb = rb - mean(rb);
    r = sum(ra .* rb) / sqrt(sum(ra.^2) * sum(rb.^2));
end
