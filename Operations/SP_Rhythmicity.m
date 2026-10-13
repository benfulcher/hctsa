function out = SP_Rhythmicity(y, m, tau, nFreq, fLow, fHigh, smoothWin, minCycles)
% SP_Rhythmicity   Phase-persistence ("rhythmicity") spectrum of a time series.
%
% A power spectrum cannot tell a sustained oscillation apart from a train of
% brief, high-power bursts: both deposit energy at the same frequency. This
% operation instead asks, at each frequency, how *predictable* the phase is a
% fixed number of cycles into the future -- high for an oscillation that keeps
% ticking (sustained/"rhythmic"), low for one that only flickers on
% intermittently (transient/"bursty").
%
% Method (the Lagged-Angle Vector Index, LAVI, of G. Karvat et al., "Universal
% rhythmic architecture uncovers two modes of neural dynamics", Nature
% Communications 17 (2026); a phase-only refinement of the lagged coherence of
% J. Fransen et al., NeuroImage 111, 158 (2015)):
%
% (1) Complex Morlet wavelet transform at each analysis frequency f (m cycles
%     per wavelet), giving the instantaneous phase phi(t,f). No Wavelet
%     Toolbox dependency -- the transform is an FFT-domain convolution here.
% (2) LAVI(f) = | mean_t exp( 1i*[ phi(t,f) - phi(t+L,f) ] ) |, with the lag
%     L = round(tau/f) samples a fixed tau cycles: the mean resultant length
%     of the lagged phase differences, in [0,1]. 1 => the phase advances by a
%     perfectly repeatable amount every tau cycles (a stable rhythm).
% (3) The spectral null. Part of LAVI is set by the power spectrum alone: a
%     Gaussian process with the same spectrum (equivalently, the mean over
%     phase-randomized surrogates) has a lag-L correlation rho(f) of its
%     wavelet output, computed here exactly from the series' own periodogram,
%     and its expected LAVI is then the mean resultant length of the phase
%     difference of a circular complex Gaussian pair,
%         LAVI0(f) = (pi/4) rho 2F1(1/2,1/2;2;rho^2)
%                  = [E(rho^2) - (1 - rho^2) K(rho^2)] / rho,
%     with K and E the complete elliptic integrals. This includes a floor that
%     the method has even for white noise: the wavelet overlaps itself across
%     the lag, so rho = exp(-pi^2 tau^2/m^2) (0.41 for the defaults) and LAVI0 =
%     0.33. The excess, LAVI(f) - LAVI0(f), is the part of the rhythmicity
%     spectrum that the power spectrum does not determine: about zero for any
%     linear Gaussian process, positive for a sustained oscillation, negative
%     for bursts. This takes the place of the surrogate "noise ribbon" of
%     Karvat et al., deterministically and at the cost of one extra FFT per
%     frequency.
%
% (The amplitude-weighted lagged coherence, |sum_t x(t,f) conj(x(t+L,f))| /
% sqrt(sum|x(t,f)|^2 sum|x(t+L,f)|^2), is deliberately not reported: up to edge
% effects it is the normalized lag-L autocorrelation of the wavelet output, a
% fixed function of the periodogram -- exactly rho(f) above.)
%
% Following hctsa convention the series is treated as evenly sampled with unit
% spacing: frequencies are in cycles per sample (0 to 0.5) and reported
% angular frequencies are 2*pi*f, as in SP_Summaries. LAVI is biased upward
% when few cycles are available, so a frequency is used only when the samples
% free of edge effects span at least minCycles cycles; the operation returns a
% single NaN if fewer than 10 frequencies qualify (for the defaults, N below
% about 240).
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
% m, width of the Morlet wavelet in cycles (default: 5, as in Karvat et al.).
%    Larger m = sharper frequency resolution, blunter time resolution.
%
% tau, the lag at which phase persistence is measured, in cycles (default:
%      1.5, as in Karvat et al., who report that the relative shape of the
%      rhythmicity spectrum is nearly independent of this choice).
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
%       (default: 5; 1 disables it). Only the band count (lavi_numPeaks) is
%       taken from the smoothed spectrum; every other statistic uses the raw
%       spectrum.
%
% minCycles, the minimum number of cycles of f that the edge-free samples must
%       span for f to be used (default: 40, which keeps the upward bias of
%       LAVI to about 0.02).
%
% ---OUTPUTS: computed over the frequencies that qualify.
% The rhythmicity spectrum LAVI(f):
% lavi_mean, lavi_median, lavi_std, lavi_iqr, lavi_max, lavi_min, lavi_range,
%       basic level and spread of rhythmicity across frequency (lavi_median is
%       the baseline).
% lavi_maxOnMedian, peak rhythmicity relative to the baseline.
% lavi_freqAtMax, lavi_freqAtMin, angular frequency (2*pi*f) of the most
%       sustained and most transient band.
% lavi_meanAboveMedian, mean amount by which above-baseline frequencies exceed
%       the baseline.
% lavi_fracSustained, lavi_fracTransient, fraction of frequencies more than 2
%       scaled-MAD above / below the baseline.
% lavi_areaSustained, area of the spectrum above its baseline, integrated over
%       log-frequency and normalized by the log-frequency span.
% lavi_ac1, lag-1 autocorrelation of the spectrum across frequency.
% lavi_numPeaks, number of local maxima of the smoothed spectrum that exceed
%       the baseline by at least one scaled-MAD (distinct rhythmic bands).
% lavi_atPeakPower, LAVI at the frequency with the most wavelet power.
% The excess over the spectral null, LAVI(f) - LAVI0(f):
% ex_mean, ex_std, ex_max, ex_min, level and spread of the rhythmicity that the
%       power spectrum does not account for.
% ex_freqAtMax, ex_freqAtMin, angular frequency of the largest and smallest
%       excess.
% ex_atPeakPower, the excess at the frequency with the most wavelet power: is
%       the dominant oscillation more sustained (> 0) or more bursty (< 0) than
%       its spectrum implies?

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
if nargin < 8 || isempty(minCycles)
    minCycles = 40;
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
if max(y) == min(y)
    out = NaN; return % constant series carries no phase information
end
y = y - mean(y);

freqs = logspace(log10(fLow), log10(fHigh), nFreq)';

% ------------------------------------------------------------------------------
%% At each frequency: wavelet transform (FFT-domain convolution; no toolbox),
%% LAVI, and its spectral null
% ------------------------------------------------------------------------------
lavi = nan(nFreq, 1);
laviNull = nan(nFreq, 1);
power = nan(nFreq, 1);
for fi = 1:nFreq
    f = freqs(fi);
    sigmaT = m / (2*pi*f); % Gaussian envelope s.d. in samples
    h = ceil(3.5 * sigmaT); % wavelet support: +/- 3.5 s.d.
    L = max(1, round(tau / f)); % lag in samples
    lo = h + 1; % samples lo:hi (and lo+L:hi+L) are free of edge effects
    hi = N - h - L;
    if (hi - lo + 1) * f < minCycles
        continue % too few cycles for an unbiased estimate: leave this frequency NaN
    end

    t = (-h:h)';
    psi = exp(2i*pi*f*t) .* exp(-t.^2 / (2*sigmaT^2));
    psi = psi - mean(psi); % zero-mean (wavelet admissibility)
    psi = psi / sqrt(sum(abs(psi).^2)); % unit energy

    % Linear convolution via FFT, then take the length-N 'same' segment:
    Lconv = N + numel(psi) - 1;
    Y = fft(y, Lconv);
    W = fft(psi, Lconv);
    c = ifft(Y .* W);
    X = c(h+1 : h+N);

    x1 = X(lo:hi);
    x2 = X(lo+L:hi+L);
    lavi(fi) = abs(mean(exp(1i * (angle(x1) - angle(x2)))));
    power(fi) = mean(abs(x1).^2);

    % Spectral null: the lag-L autocorrelation of the wavelet output implied by
    % the (zero-padded, so linear) periodogram, and the LAVI it implies for a
    % Gaussian process:
    r = ifft(abs(Y).^2 .* abs(W).^2);
    laviNull(fi) = SUB_laviFromRho(abs(r(L+1)) / real(r(1)));
end

% ------------------------------------------------------------------------------
%% Summarize
% ------------------------------------------------------------------------------
valid = isfinite(lavi) & isfinite(laviNull);
if sum(valid) < 10
    out = NaN; return % too little of the spectrum is estimable to summarize
end

f = freqs(valid);
logf = log10(f);
S = lavi(valid);
ex = S - laviNull(valid);
[~, iPow] = max(power(valid));

% The rhythmicity spectrum LAVI(f):
med = median(S);
scaledMAD = 1.4826 * median(abs(S - med));

out.lavi_mean = mean(S);
out.lavi_median = med;
out.lavi_std = std(S);
out.lavi_iqr = iqr(S);
out.lavi_max = max(S);
out.lavi_min = min(S);
out.lavi_range = max(S) - min(S);
out.lavi_maxOnMedian = max(S) / med;

[~, iMax] = max(S);
[~, iMin] = min(S);
out.lavi_freqAtMax = 2*pi*f(iMax);
out.lavi_freqAtMin = 2*pi*f(iMin);

above = S(S > med) - med;
if isempty(above)
    out.lavi_meanAboveMedian = 0;
else
    out.lavi_meanAboveMedian = mean(above);
end

if scaledMAD < eps
    out.lavi_fracSustained = 0;
    out.lavi_fracTransient = 0;
else
    out.lavi_fracSustained = mean(S > med + 2*scaledMAD);
    out.lavi_fracTransient = mean(S < med - 2*scaledMAD);
end

% Area above baseline, integrated over log-frequency, span-normalized:
out.lavi_areaSustained = trapz(logf, max(S - med, 0)) / (logf(end) - logf(1));

% Lag-1 autocorrelation of the spectrum across frequency:
S0 = S - mean(S);
denom = sum(S0.^2);
if denom < eps
    out.lavi_ac1 = NaN;
else
    out.lavi_ac1 = sum(S0(1:end-1) .* S0(2:end)) / denom;
end

% Number of rhythmic bands: local maxima of the smoothed spectrum that clear
% the baseline by at least one scaled-MAD (a Daniell/moving-average window
% stands in for Welch-style segment averaging, which the series is too short
% to afford here).
Ssm = movmean(S, smoothWin);
dS = diff(Ssm);
isPeak = [false; (dS(1:end-1) > 0) & (dS(2:end) <= 0); false];
out.lavi_numPeaks = sum(isPeak & (Ssm > med + scaledMAD));

out.lavi_atPeakPower = S(iPow);

% The excess over the spectral null:
out.ex_mean = mean(ex);
out.ex_std = std(ex);
out.ex_max = max(ex);
out.ex_min = min(ex);
[~, iMax] = max(ex);
[~, iMin] = min(ex);
out.ex_freqAtMax = 2*pi*f(iMax);
out.ex_freqAtMin = 2*pi*f(iMin);
out.ex_atPeakPower = ex(iPow);

end

% ==============================================================================
function g = SUB_laviFromRho(rho)
    % Mean resultant length of the phase difference of a circular complex
    % Gaussian pair with correlation magnitude rho:
    % (pi/4) rho 2F1(1/2,1/2;2;rho^2) = [E(rho^2) - (1 - rho^2) K(rho^2)] / rho.
    if rho < 1e-6
        g = pi * rho / 4; % small-rho limit (avoids 0/0)
    elseif rho >= 1
        g = 1;
    else
        [K, E] = ellipke(rho^2);
        g = (E - (1 - rho^2) * K) / rho;
    end
end
