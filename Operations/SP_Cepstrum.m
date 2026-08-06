function out = SP_Cepstrum(y, maxPeriod, minPeriod)
% SP_Cepstrum   Cepstral statistics: harmonic (comb) structure of the power spectrum
%
% Computes the real cepstrum, the inverse Fourier transform of the log
% magnitude spectrum, and summarizes the structure of its dominant peak.
%
% The cepstrum answers a question none of hctsa's other spectral
% operations ask: whether the peaks in the spectrum are *harmonically
% related* to one another. SP_Summaries counts and characterizes spectral
% peaks (numPeaks, maxProm, peakPower_*), but those statistics are
% identical for a signal with a fundamental plus five harmonics and one
% with five arbitrarily-placed peaks. A harmonic series is a comb of
% peaks evenly spaced by the fundamental frequency f0, which is itself a
% periodicity *in the spectrum*; taking the spectrum of the log spectrum
% collapses that whole comb into a single cepstral peak at a quefrency
% ('frequency' spelled backwards -- the cepstral axis has units of time)
% equal to the fundamental period 1/f0.
%
% Note SP_Summaries' spectral autocorrelation fields (ac1/ac2/tau) do not
% cover this: they are evaluated at lags of 1-2 frequency bins and at the
% ACF's first zero crossing, all far too short-lag to detect harmonic
% spacing, which is typically tens of bins. The cepstrum is also a
% different proposition from time-domain periodicity operations
% (PD_PeriodicityWang, CO_AutoCorr on the raw series): it is computed
% from the log spectrum, so it is blind to phase and responds to
% multiplicative rather than additive structure, which is why it remains
% sensitive to a harmonic series sitting on top of a strong 1/f
% background.
%
% ---INPUTS:
% y, the input time series
%
% maxPeriod, the longest fundamental period (in samples) to search for
%            (default: 100). Deliberately a fixed number of samples
%            rather than a fraction of the series length: making the
%            search range scale with N would make the resulting
%            quefrency statistics partly a restatement of N rather than
%            of the data (cf. the N-dependence bug found in
%            SP_Summaries' tau). Series shorter than 4*maxPeriod (i.e.
%            too short to contain several cycles of the longest period
%            searched) return NaN rather than an unreliable estimate.
%
% minPeriod, the shortest fundamental period (in samples) to search for
%            (default: 4). The low-quefrency end of the cepstrum encodes
%            the smooth spectral envelope (overall tilt / 1-f slope)
%            rather than harmonic structure, so it is excluded. A period
%            of 4 samples is also about the shortest that can support a
%            harmonic series at all: the second harmonic of f0 = 0.25
%            cycles/sample already sits at the Nyquist frequency.
%
% ---WHAT THIS DOES AND DOESN'T DETECT (validated on synthetic signals):
% Because it detects a *comb*, this operation responds to harmonic
% richness, not to periodicity as such -- which is what makes it
% complementary to time-domain periodicity operations rather than a
% duplicate of them:
%   - A pure sinusoid, however strongly periodic, has no harmonics and
%     hence no comb: it scores at the noise floor here (CPP ~ 0.05
%     against a white-noise null of ~0.025), while a time-domain
%     autocorrelation would flag it very strongly.
%   - Detecting a fundamental requires enough harmonics spread across
%     the band. A period of 64 samples with only 5 harmonics (which
%     span just the lowest 8% of the spectrum) is missed; the same
%     period with 10 or more harmonics is recovered exactly. CPP rises
%     monotonically with the number of comb teeth.
%   - Failure is quiet rather than confidently wrong: in the missed case
%     above, CPP stays near the null (0.064, against 0.10-0.18 for
%     successfully-detected combs), so a weak CPP correctly indicates
%     'no comb found' rather than endorsing a spurious period. The
%     `period` output is only meaningful when CPP/peakRatio indicate a
%     genuine peak.
%
% ---OUTPUTS:
% The estimated fundamental period (quefrency of the dominant cepstral
% peak) and the height of that peak; the cepstral peak prominence (CPP),
% the standard robust measure, being the peak height above a linear
% regression fit through the cepstrum across the search range (this
% normalizes away the overall cepstral trend, so it does not simply
% track the spectrum's dynamic range); the peak height in units of the
% standard deviation of the cepstrum over the search range; the
% rahmonic ratio, comparing the cepstrum at twice the peak quefrency to
% the peak itself (a genuine harmonic comb repeats at multiples of the
% fundamental period, so a real periodicity shows a secondary 'rahmonic'
% peak, whereas an isolated fluke does not); and the mean and standard
% deviation of the cepstrum over the search range.
%
% Of these, period/CPP/peakRatio/rahmonicRatio/meanCeps are registered as
% hctsa features. `peak` and `stdCeps` are computed (CPP and peakRatio
% are derived from them) but not registered: checked across 1000 real
% series, raw `peak` is a near-duplicate of CPP (R^2 = 0.98) and the
% baseline-subtracted CPP is both the standard published formulation and
% the better-behaved one, since it does not simply track the cepstrum's
% overall scale; `stdCeps` is likewise largely redundant with CPP
% (R^2 = 0.74) and is the most ad hoc quantity of the set, its real role
% being as the normalizer inside peakRatio.
%
% On novelty: across the same 1000 series these fields reach a maximum
% R^2 of only 0.52 against all ~7700 other hctsa features, and 0.20
% against the ~690 periodicity/autocorrelation-adjacent features
% specifically -- i.e. harmonic-comb structure really was not being
% captured elsewhere in the library. rahmonicRatio is the most
% independent single feature found (max R^2 0.04 against anything).

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

if nargin < 2 || isempty(maxPeriod)
    maxPeriod = 100;
end
if nargin < 3 || isempty(minPeriod)
    minPeriod = 4;
end

if minPeriod < 2
    error('minPeriod = %u is below the Nyquist limit (a period needs >= 2 samples)', minPeriod);
end
if maxPeriod <= minPeriod
    error('maxPeriod (%u) must exceed minPeriod (%u)', maxPeriod, minPeriod);
end

N = length(y);
minCycles = 4; % need several cycles of the longest period searched for a meaningful estimate
if N < minCycles * maxPeriod
    warning(['Time series (N = %u) too short to search for periods up to %u samples ' ...
             '(need >= %u)'], N, maxPeriod, minCycles * maxPeriod);
    out = NaN; return
end

if all(y == y(1)) % constant series has an all-zero spectrum -> log(0)
    warning('Constant time series has no spectral (or cepstral) structure');
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Real cepstrum
% ------------------------------------------------------------------------------
% c(q) = IFFT{ log|FFT(y)| }. Taking the log *before* the second
% transform is what makes this multiplicative-structure-sensitive: a
% harmonic comb multiplies the spectral envelope, and the log turns that
% product into a sum, so the comb separates from the envelope into a
% distinct (high-quefrency) part of the cepstrum.
NFFT = 2^nextpow2(N);
X = fft(y, NFFT);
logMag = log(abs(X) + eps); % eps guards spectral nulls (|X| exactly 0)

% Remove the smooth spectral envelope before transforming. Without this,
% the low-quefrency end of the cepstrum is dominated by the overall
% shape of the spectrum (its tilt / 1-f slope) rather than by harmonic
% structure, and for strongly-coloured spectra that envelope leaks far
% enough up the quefrency axis to produce outright false positives:
% before this step, pure 1/f noise -- which has no periodicity
% whatsoever -- scored a *higher* peak-to-spread ratio (5.35) than a
% genuine 5-harmonic comb (3.60), and a true period of 64 samples was
% missed entirely in favour of an envelope artefact at low quefrency.
% A low-order polynomial in frequency is exactly the 'smooth envelope'
% component, so subtracting a fitted one leaves the comb ripple that the
% cepstrum is meant to detect. Order 4 is enough to absorb a 1-f-style
% tilt and gentle curvature without being flexible enough to start
% fitting the ripple itself.
envOrder = 4;
nHalf = floor(NFFT / 2) + 1;
halfLogMag = logMag(1:nHalf);
fIdx = (0:nHalf - 1)' / (nHalf - 1); % normalized frequency axis for conditioning
pEnv = polyfit(fIdx, halfLogMag, envOrder);
halfDetrended = halfLogMag - polyval(pEnv, fIdx);

% Mirror back to a full Hermitian-symmetric spectrum so the cepstrum is real:
logMagDetrended = [halfDetrended; flipud(halfDetrended(2:end - 1))];
c = real(ifft(logMagDetrended));

% Quefrency index q (1-based) corresponds to a period of (q-1) samples,
% so a period of p samples sits at index p+1:
qIdx = (minPeriod:maxPeriod)' + 1;
if qIdx(end) > floor(NFFT / 2)
    % Shouldn't be reachable given the length check above, but the
    % cepstrum is only meaningful over its first half (it is symmetric):
    qIdx = qIdx(qIdx <= floor(NFFT / 2));
end
cSearch = c(qIdx);
periods = qIdx - 1; % period in samples corresponding to each search index

% ------------------------------------------------------------------------------
%% Dominant cepstral peak
% ------------------------------------------------------------------------------
[peakVal, iPeak] = max(cSearch);
out.period = periods(iPeak); % estimated fundamental period, in samples
out.peak = peakVal;

% Basic distributional context over the search range:
out.meanCeps = mean(cSearch);
out.stdCeps = std(cSearch);

% Peak height in units of the cepstrum's own spread over the search range
% (scale-free, unlike `peak` itself):
if out.stdCeps > 0
    out.peakRatio = (peakVal - out.meanCeps) / out.stdCeps;
else
    out.peakRatio = NaN;
end

% ------------------------------------------------------------------------------
%% Cepstral peak prominence (CPP)
% ------------------------------------------------------------------------------
% The standard robust formulation: the peak's height above a linear
% regression through the cepstrum over the search range. Subtracting the
% fitted trend removes the overall decay of the cepstrum with quefrency,
% so CPP measures how much the peak stands out from its local
% background rather than how large the cepstrum is overall.
pFit = polyfit(periods, cSearch, 1);
baseline = polyval(pFit, periods);
out.CPP = peakVal - baseline(iPeak);

% ------------------------------------------------------------------------------
%% Rahmonic structure
% ------------------------------------------------------------------------------
% A true fundamental period p produces cepstral energy at 2p, 3p, ...
% ('rahmonics'), because a harmonic comb in the spectrum is itself
% periodic. An isolated spurious peak does not. Compare the cepstrum at
% twice the peak quefrency to the peak, both measured against the fitted
% baseline so the comparison isn't dominated by the cepstral trend:
q2 = 2 * out.period + 1;
peakAboveBase = peakVal - baseline(iPeak);
if q2 <= floor(NFFT / 2) && peakAboveBase > 0
    rahmonicAboveBase = c(q2) - polyval(pFit, 2 * out.period);
    out.rahmonicRatio = rahmonicAboveBase / peakAboveBase;
else
    out.rahmonicRatio = NaN;
end

end
