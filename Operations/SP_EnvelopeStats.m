function out = SP_EnvelopeStats(y, powerFrac, trimFrac)
% SP_EnvelopeStats   Statistics of the amplitude envelope of the full series and of its dominant oscillation.
%
% Computes the instantaneous amplitude envelope |z(t)| of the analytic signal
% z(t) = y(t) + i H[y](t), where H is the Hilbert transform, for (a) the full
% series and (b) the dominant band: the narrowest frequency band centered on the
% largest periodogram peak (DC and Nyquist bins excluded) that holds a fraction
% powerFrac of the total power, and at least 2 bins either side of the peak. Its
% width therefore adapts to the series: it is the width of the dominant peak
% for a narrowband oscillation, and a large part of the spectrum for broadband noise.
% Each analytic signal comes from the FFT: zeroing the bins outside the band
% (negative frequencies included) and doubling the rest before an inverse FFT,
% with no toolbox needed.
%
% The envelope summaries describe amplitude modulation: how variable the
% envelope is (coefficient of variation), how asymmetric and heavy-tailed its
% distribution is (kurtosis for the full series; skewness and kurtosis for the
% dominant band), and how long it takes to decorrelate (the 1/e timescale of its
% autocorrelation function). An unmodulated sinusoid has a constant envelope (CV
% near 0); bursting or amplitude-modulated signals have a large CV, and a high
% kurtosis for intermittent bursts.
%
% Baseline for Gaussian noise: the analytic signal of a stationary Gaussian
% process is complex Gaussian, so its envelope is Rayleigh distributed:
% CV = sqrt(4/pi - 1) = 0.5227, skewness 0.6311, kurtosis 3.2451 (checked
% numerically on white noise for the full-band fields; the full-band skewness is
% not returned, as it is redundant with the kurtosis). Values of the CV below
% this indicate an envelope steadier than noise (e.g., a sinusoid in noise follows
% a Rice distribution), and above it an envelope more modulated than noise.
%
% To limit edge effects (the FFT filter is circular, so the series ends wrap
% around), trimFrac of the samples are dropped from each end of the envelope and
% phase before any summary is computed. Timescales are in samples (as elsewhere in
% hctsa), so they scale with the sampling rate (for a narrowband oscillation, the
% dominant-band timescale and frequency spread do too, because the band is set by
% the width of the spectral peak). The dominant-band fields of a narrow band are
% based on few effectively independent envelope values (about N times the band
% width), so they are biased toward lower CV for short series.
%
% ---INPUTS:
% y, the input time series
% powerFrac, the fraction of the total (one-sided, DC and Nyquist excluded) spectral power
%            that the dominant band, centered on the largest periodogram peak, must
%            contain (default: 0.5)
% trimFrac, the fraction of samples dropped from each end of the analytic signal
%           before computing summaries (default: 0.05)
%
% ---OUTPUTS:
% full_cv, the envelope's coefficient of variation (standard deviation over mean), full series
% full_kurt, the envelope's kurtosis, full series
% full_tau, the 1/e decay time (in samples) of the envelope's autocorrelation function, full series
% dom_cv, dom_skew (the envelope's skewness), dom_kurt, dom_tau, as above for
%         the dominant band
% dom_ifspread, a robust spread of the dominant band's instantaneous frequency
%               (1.4826 times the median absolute deviation of the phase
%               increments, in cycles per sample)
% All fields are NaN for constant, non-finite, or very short (N < 50) series.
% A 1/e timescale is NaN when the autocorrelation never falls below 1/e within N/2 lags.
% For an envelope that is essentially constant (as for a sinusoid without noise)
% the skewness, kurtosis and timescale are an undefined 0/0 (they divide by the
% variance), and their limit depends on how the envelope becomes constant (a
% sinusoid with vanishing added noise and one with vanishing amplitude modulation
% tend to different values). They are set by convention to skewness 0 and
% kurtosis 3 (the Gaussian values) and a 1/e timescale equal to the largest lag
% searched, floor(n/2) for the n samples left after trimming (a constant envelope
% never decorrelates). So that no threshold is needed, each measured value is
% shrunk toward its conventional value with the weight w = c^2/(c^2 + 1e-20),
% where c is the envelope's CV (computed as usual, and reported unchanged):
% reported = convention + w*(measured - convention). The reported values vary
% continuously with c; the change happens over c of about 1e-11 to 1e-9, and they
% are unchanged (to a relative 1e-12 or better) when c exceeds 1e-4. The scale
% 1e-10 is four to five orders above the round-off of a noiseless sinusoid (c of
% order 1e-15 to 1e-14) and far below any real variability. The weights are
% applied to the full and dominant-band envelopes separately.
%
% ---REFERENCES:
% B. Boashash, "Estimating and interpreting the instantaneous frequency of a
% signal. I. Fundamentals", Proc. IEEE 80(4), 520-538 (1992).
%
% ---NOTES:
% The dominant band is the narrowest window around the largest periodogram peak
% holding half the power (by default), so for a narrowband oscillation it follows the
% width of the spectral peak and dom_tau and dom_ifspread scale with the
% sampling rate like full_tau (decimating by 2 and 4 gave dom_tau ratios of
% about 1/2 and 1/4 and dom_ifspread ratios of 2 and 4, in simulation). A peak
% narrower than the frequency resolution (a sinusoid in noise, a random walk)
% gives a band of the 2-bin minimum, so dom_tau then grows in proportion to N. For
% broadband noise the band is a large part of the spectrum, so dom_tau is only a few
% samples.
%
% Also computed during development but not kept, as redundant: the lag-1
% envelope autocorrelation (r = 0.95 with the 1/e timescale for the full band;
% near 1 for every series for a narrow band), and the instantaneous-frequency
% spread relative to its median (r = 0.92 with SP_PhaseFluctuationScaling's meanFreq).

% ------------------------------------------------------------------------------
% Copyright (C) 2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
%% Check inputs, set defaults
% ------------------------------------------------------------------------------
y = y(:);
if nargin < 2 || isempty(powerFrac)
    powerFrac = 0.5;
end
if nargin < 3 || isempty(trimFrac)
    trimFrac = 0.05;
end

% Fixed set of output fields (NaN whenever a value is undefined):
fldFull = {'full_cv', 'full_kurt', 'full_tau'};
fldDom = {'dom_cv', 'dom_skew', 'dom_kurt', 'dom_tau', 'dom_ifspread'};
for f = [fldFull, fldDom]
    out.(f{1}) = NaN;
end

N = length(y);
if N < 50 || any(~isfinite(y)) || std(y) == 0
    return
end
y = y - mean(y);

% ------------------------------------------------------------------------------
%% Frequency bins (DC and Nyquist excluded, as in SP_PhaseAmpCoupling)
% ------------------------------------------------------------------------------
halfN = floor(N / 2) + 1;
if mod(N, 2) == 0
    usableBins = 2:(halfN - 1);
else
    usableBins = 2:halfN;
end
Y = fft(y);

% Samples dropped from each end of the (circularly computed) analytic signal:
nTrim = max(1, round(trimFrac * N));
keep = (nTrim + 1):(N - nTrim);

% ------------------------------------------------------------------------------
%% (a) Full band
% ------------------------------------------------------------------------------
Yfull = zeros(N, 1);
Yfull(usableBins) = 2 * Y(usableBins);
zFull = ifft(Yfull);
envFull = abs(zFull(keep));
[out.full_cv, ~, out.full_kurt, out.full_tau] = envelopeSummary(envFull);

% ------------------------------------------------------------------------------
%% (b) Dominant band: the largest periodogram peak +/- the smallest half-width
% (at least 2 bins) for which the band holds a fraction powerFrac of the total power
% ------------------------------------------------------------------------------
pow = abs(Y(usableBins)).^2;
nBins = length(usableBins);
[~, iPeak] = max(pow);
cumPow = [0; cumsum(pow)];
hws = (2:nBins)';
bandPow = cumPow(min(nBins, iPeak + hws) + 1) - cumPow(max(1, iPeak - hws));
halfWidthBins = hws(find(bandPow >= powerFrac * cumPow(end), 1));
if isempty(halfWidthBins)
    halfWidthBins = nBins;
end
peakBin = usableBins(iPeak);
bandBins = max(usableBins(1), peakBin - halfWidthBins):min(usableBins(end), peakBin + halfWidthBins);

Ydom = zeros(N, 1);
Ydom(bandBins) = 2 * Y(bandBins);
zDom = ifft(Ydom);
envDom = abs(zDom(keep));
[out.dom_cv, out.dom_skew, out.dom_kurt, out.dom_tau] = envelopeSummary(envDom);

% Instantaneous frequency (cycles per sample): the unwrapped phase increments
% over the trimmed segment, summarized robustly (scaled MAD, which equals the
% standard deviation for a Gaussian)
phi = unwrap(angle(zDom(keep)));
instFreq = diff(phi) / (2 * pi);
out.dom_ifspread = 1.4826 * median(abs(instFreq - median(instFreq)));

end

% ------------------------------------------------------------------------------
function [cv, sk, ku, tau] = envelopeSummary(env)
% Distributional and autocorrelation summaries of an amplitude envelope
cv = NaN; sk = NaN; ku = NaN; tau = NaN;
n = length(env);
m = mean(env);
if ~(m > 0)
    return
end
e = env - m;
s2 = mean(e.^2);
cv = sqrt(s2) / m; % population standard deviation over the mean

% Weight of the measured skewness, kurtosis and timescale: these are a 0/0 for a
% constant envelope, shrunk continuously toward conventional values (see the
% function help); w -> 1 for any envelope that varies
cvScale = 1e-10; % well above the round-off CV of a noiseless sinusoid (~1e-15)
w = cv^2 / (cv^2 + cvScale^2);
tauMax = floor(n / 2); % the largest lag searched below (the envelope never decorrelates)

skMeas = NaN; kuMeas = NaN; tauMeas = NaN;
if s2 > 0
    skMeas = mean(e.^3) / s2^1.5;
    kuMeas = mean(e.^4) / s2^2;

    % Autocorrelation of the envelope via the FFT (zero-padded, biased estimator)
    nfft = 2^nextpow2(2 * n);
    F = fft(e, nfft);
    acf = real(ifft(abs(F).^2));
    acf = acf(1:floor(n / 2) + 1) / acf(1); % lags 0..N/2
    iCross = find(acf < exp(-1), 1);
    if ~isempty(iCross) && iCross > 1
        % linear interpolation between lags (iCross-2) and (iCross-1)
        a0 = acf(iCross - 1);
        a1 = acf(iCross);
        tauMeas = (iCross - 2) + (a0 - exp(-1)) / (a0 - a1);
    end
end

sk = shrinkToLimit(skMeas, 0, w);
ku = shrinkToLimit(kuMeas, 3, w);
tau = shrinkToLimit(tauMeas, tauMax, w);
end

% ------------------------------------------------------------------------------
function v = shrinkToLimit(measured, limit, w)
% limit + w*(measured - limit), where limit is the conventional value; an
% undefined measured value takes it only when the weight w is negligible (an
% essentially constant envelope)
if isnan(measured)
    if w < 1e-6
        v = limit;
    else
        v = NaN;
    end
else
    v = limit + w * (measured - limit);
end
end
