function out = SP_Rhythmicity(y, m, tau, nFreq, fLow, fHigh)
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
%     per wavelet), giving an instantaneous phase phi(t,f). No Wavelet Toolbox
%     dependency -- the transform is an FFT-domain convolution here.
% (2) At each f, take the phase difference between phi(t,f) and phi(t+L,f),
%     where the lag L = round(tau/f) samples is a fixed tau cycles. LAVI(f) is
%     the mean resultant length of that phase-difference distribution:
%       LAVI(f) = | mean_t exp( 1i*[ phi(t,f) - phi(t+L,f) ] ) |,   in [0,1].
%     1 => the phase advances by a perfectly repeatable amount every tau
%     cycles (a stable rhythm); ~0 => the phase L samples ahead is
%     unpredictable (noise, or a rhythm too short-lived to bridge the lag).
% (3) The rhythmicity spectrum LAVI(f) is summarized below. Unlike a power
%     spectrum it has a stable baseline (its own median across frequency),
%     so peaks above baseline flag genuinely sustained bands and troughs
%     below it flag transient/bursty bands.
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
% ---OUTPUTS: statistics of the rhythmicity spectrum LAVI(f), computed over
% the frequencies that survive the length check:
% lavi_mean, lavi_median, lavi_std, lavi_iqr, lavi_max, lavi_min, lavi_range,
%       basic level and spread of rhythmicity across frequency. lavi_median is
%       the natural baseline (Karvat et al.).
% lavi_maxOnMedian, peak rhythmicity relative to the baseline -- how far the
%       most sustained band stands out.
% lavi_freqAtMax, lavi_freqAtMin, angular frequency (2*pi*f) of the most
%       sustained and most transient band.
% lavi_meanAboveMedian, mean amount by which above-baseline frequencies exceed
%       the baseline (average prominence of sustained bands).
% lavi_fracSustained, lavi_fracTransient, fraction of frequencies sitting more
%       than 2 scaled-MAD above / below the baseline -- a data-driven stand-in
%       for the surrogate "noise ribbon" of the original method (which needs
%       many minutes of data and is not attempted here).
% lavi_areaSustained, area of LAVI(f) above its baseline, integrated over
%       log-frequency and normalized by the log-frequency span: the total
%       "sustained mass" of the spectrum.
% lavi_slope, ordinary-least-squares slope of LAVI(f) against log10(f): does
%       rhythmicity concentrate at low or high frequencies?
% lavi_ac1, lag-1 autocorrelation of the LAVI(f) sequence: a smooth profile
%       (broad bands) vs. a jagged one.
% lavi_numPeaks, number of local maxima of LAVI(f) that exceed the baseline
%       (number of distinct sustained bands).
% lavi_corrPower, Spearman correlation between LAVI(f) and log10 of the mean
%       wavelet power at f: near 0 => rhythmicity carries information beyond
%       the power spectrum; near 1 => it largely tracks power.

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
% Per frequency we also record the wavelet half-length, so the lagged-phase
% step below can discard the edge-contaminated samples.
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
%% Lagged-angle vector index (LAVI) at each frequency
% ------------------------------------------------------------------------------
lavi = nan(nFreq, 1);
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

    ph1 = angle(X(lo:hi, fi));
    ph2 = angle(X(lo+L:hi+L, fi));
    lavi(fi) = abs(mean(exp(1i * (ph1 - ph2))));
    logPower(fi) = log10(mean(abs(X(lo:hi, fi)).^2));
end

% ------------------------------------------------------------------------------
%% Summarize the rhythmicity spectrum
% ------------------------------------------------------------------------------
valid = isfinite(lavi);
if sum(valid) < 10
    out = NaN; return % too little of the spectrum is estimable to summarize
end

L = lavi(valid);
f = freqs(valid);
logf = log10(f);
P = logPower(valid);
med = median(L);

out.lavi_mean = mean(L);
out.lavi_median = med;
out.lavi_std = std(L);
out.lavi_iqr = iqr(L);
out.lavi_max = max(L);
out.lavi_min = min(L);
out.lavi_range = max(L) - min(L);
out.lavi_maxOnMedian = max(L) / med;

[~, iMax] = max(L);
[~, iMin] = min(L);
out.lavi_freqAtMax = 2*pi*f(iMax);
out.lavi_freqAtMin = 2*pi*f(iMin);

above = L(L > med) - med;
if isempty(above)
    out.lavi_meanAboveMedian = 0;
else
    out.lavi_meanAboveMedian = mean(above);
end

scaledMAD = 1.4826 * median(abs(L - med));
if scaledMAD < eps
    out.lavi_fracSustained = 0;
    out.lavi_fracTransient = 0;
else
    out.lavi_fracSustained = mean(L > med + 2*scaledMAD);
    out.lavi_fracTransient = mean(L < med - 2*scaledMAD);
end

% Area above baseline, integrated over log-frequency, span-normalized:
out.lavi_areaSustained = trapz(logf, max(L - med, 0)) / (logf(end) - logf(1));

% OLS slope of rhythmicity against log-frequency:
pf = polyfit(logf, L, 1);
out.lavi_slope = pf(1);

% Lag-1 autocorrelation of the LAVI(f) sequence:
L0 = L - mean(L);
denom = sum(L0.^2);
if denom < eps
    out.lavi_ac1 = NaN;
else
    out.lavi_ac1 = sum(L0(1:end-1) .* L0(2:end)) / denom;
end

% Number of above-baseline local maxima (distinct sustained bands):
isPeak = [false; (L(2:end-1) > L(1:end-2)) & (L(2:end-1) > L(3:end)); false];
out.lavi_numPeaks = sum(isPeak & (L > med));

% Does rhythmicity just track power?
if std(P) < eps
    out.lavi_corrPower = NaN;
else
    out.lavi_corrPower = subSpearman(L, P);
end

% ------------------------------------------------------------------------------
function r = subSpearman(a, b)
    % Spearman rank correlation without a Statistics Toolbox dependency
    % (continuous wavelet quantities: ties are measure-zero, ignored).
    n = numel(a);
    [~, oa] = sort(a); ra = zeros(n,1); ra(oa) = 1:n;
    [~, ob] = sort(b); rb = zeros(n,1); rb(ob) = 1:n;
    ra = ra - mean(ra);
    rb = rb - mean(rb);
    r = sum(ra .* rb) / sqrt(sum(ra.^2) * sum(rb.^2));
end

end
