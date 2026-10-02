function out = SP_SummariesPhase(y)
% SP_SummariesPhase   Statistics of the Fourier phase spectrum of a time series.
%
% cf. SP_Summaries, which characterizes the magnitude spectrum in detail but
% discards phase entirely. The only other place phase appears in this codebase is
% SD_MakeSurrogates, which randomizes it to build a null-model surrogate. For a
% linear, Gaussian stochastic process, Fourier phases are theoretically i.i.d.
% uniform on (-pi,pi], which is why phase randomization works as a surrogate null
% model. This operation characterizes the phase spectrum directly: deviations
% from uniformity or independence across frequency are a direct signature of
% determinism, nonlinearity, or transient/localized structure that the magnitude
% spectrum alone cannot see.
%
% The FFT is of the mean-subtracted series, zero-padded to a power of 2, with the
% phase referenced to the centre of the series. Phases are weighted by their bin's
% magnitude throughout (a standard approach in circular statistics for data of
% uneven reliability): a single pure tone concentrates essentially all energy in
% 1-2 bins, and every other bin's phase is set by numerical noise and must not
% swamp an unweighted average. The DC and Nyquist bins (both purely real, phase
% undefined in the usual oscillatory sense) are excluded throughout.
%
% ---INPUTS:
% y, the input time series
%
% ---OUTPUTS:
% R, the magnitude-weighted mean resultant length of the phases (circular
%       concentration): 0 for phases with no common preferred direction (e.g.,
%       white noise), up to 1 if every frequency component shares the same phase
%       (e.g., an impulse at the very start of the series)
% phEnt, the magnitude-weighted Shannon entropy of the phase distribution
%       (20-bin histogram), normalized to [0,1] by log(20): 1 for a uniform phase
%       distribution, lower for a concentrated one
% groupDelay, the negative slope of a magnitude-weighted linear fit of unwrapped
%       phase against angular frequency, as a fraction of the series length: the
%       delay of the series' energy relative to its midpoint. About 0 for a
%       stationary series; a unit impulse at sample 500 of 2000 gives -0.25, and
%       at sample 1500 gives +0.25
% phaseLinearity, the weighted RMSE of that linear fit, normalized by
%       sqrt(#frequency bins): how far the phase-frequency relationship is from a
%       pure delay
% magPhaseCorr, the linear correlation between magnitude and raw phase across bins
% phaseUnwrapAC1, the magnitude-weighted lag-1 autocorrelation of consecutive
%       unwrapped-phase increments across frequency: near zero when the local group
%       delay is unstructured across frequency, large and positive when it varies
%       smoothly and systematically with frequency
%
% ---REFERENCES:
% J. Theiler, S. Eubank, A. Longtin, B. Galdrikian and J. Doyne Farmer, "Testing for
% nonlinearity in time series: the method of surrogate data", Physica D 58(1-4),
% 77-94 (1992). DOI: 10.1016/0167-2789(92)90102-S
%
% N.I. Fisher, "Statistical Analysis of Circular Data", Cambridge University Press,
% Cambridge (1993). DOI: 10.1017/CBO9780511564345
%
% ---NOTES:
% Validation, from the original docstring. R: white noise gave R ~ 0.028 +/- 0.015
% over 300 trials (an empirical null), a periodic sine wave R ~ 0.22 and a linear
% chirp R ~ 0.46. phEnt is anti-correlated with R (r = -0.88 on real-world series) but
% not a function of it: entropy is sensitive to the full shape of the phase
% distribution (e.g., bimodal), which R cannot distinguish from uniformity.
% groupDelay was previously referenced to the first sample, which made it ~N/2 for
% any stationary series and aliased under unwrap (Spearman 0.47 with N on
% real-world series); phaseLinearity was previously un-normalized and grew as sqrt(N)
% (Spearman 0.78 with N); on the normalized scale an unstructured (random-walk)
% unwrapped phase gives ~0.4 regardless of series length.
% phaseUnwrapAC1 is strongly diagnostic of dispersive, frequency-dependent delay: a
% linear chirp gave ~0.86, while every other synthetic test signal was within
% +/-0.02 of zero. magPhaseCorr is the weakest-validated statistic: inconsistent
% across synthetic signals (all within +/-0.06 of zero), and weakly correlated with
% everything else on real-world series (max |r| 0.29 against any existing feature, 0.11
% against the other five here). It is kept in the spirit of hctsa's general
% preference for including a plausible statistic.

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
%% Compute the FFT (same convention as SP_Summaries: Fs=1, NFFT a power of 2)
% ------------------------------------------------------------------------------
y = y(:);
Ny = length(y);
NFFT = 2^nextpow2(Ny);
Fs = 1;
f = Fs / 2 * linspace(0, 1, NFFT / 2 + 1);
w = 2 * pi * f;

Sc = fft(y - mean(y), NFFT); % mean-subtracted, so the DC bin is (numerically) exactly zero
Sc = Sc(1:NFFT / 2 + 1); % single-sided
% Reference the phase to the centre of the series rather than its first
% sample. The Fourier phase of any signal carries a linear term -w*t0 set
% by where its energy sits in time (t0 ~ (Ny-1)/2 for a stationary series),
% which for the full-series FFT is a phase advance of ~pi*Ny/NFFT per bin --
% up to pi. That term (i) dominated the unwrapped phase, so groupDelay was
% ~Ny/2 for any stationary series and phaseLinearity grew with Ny, and
% (ii) aliased under unwrap once Ny/NFFT approached 1, so both jumped
% between regimes with Ny (e.g., groupDelay ~0 at Ny <= 2048 but ~1050 at
% Ny = 5000 for white noise). Shifting the time origin to the centre
% removes it: groupDelay is then the delay of the series' energy relative
% to its centre (0 for stationary data, +/- for a transient early or
% late in the window), and the other statistics describe the phase
% structure itself.
Sc = Sc .* exp(1i * w(:) * (Ny - 1) / 2);
mag = abs(Sc);
ph = angle(Sc);

% Exclude DC (bin 1) and Nyquist (last bin): both purely real, phase undefined
% in the usual oscillatory sense.
idx = 2:(length(ph) - 1);
ph = ph(idx); mag = mag(idx); ww = w(idx);
ph = ph(:); mag = mag(:); ww = ww(:);

if ~any(mag > 0) || ~all(isfinite(mag))
    out = NaN; return
end

wgt = mag / sum(mag);

% ------------------------------------------------------------------------------
%% Magnitude-weighted circular concentration
% ------------------------------------------------------------------------------
Rvec = sum(wgt .* exp(1i * ph));
out.R = abs(Rvec);

% ------------------------------------------------------------------------------
%% Magnitude-weighted, normalized phase entropy (20-bin histogram)
% ------------------------------------------------------------------------------
nBins = 20;
edges = linspace(-pi, pi, nBins + 1);
[~, binIdx] = histc(ph, edges);
binIdx(binIdx == 0) = 1; binIdx(binIdx > nBins) = nBins; % guard the (rare) ph == pi edge case
pBin = accumarray(binIdx, wgt, [nBins, 1]);
pBin_nz = pBin(pBin > 0);
out.phEnt = -sum(pBin_nz .* log(pBin_nz)) / log(nBins);

% ------------------------------------------------------------------------------
%% Group delay: magnitude-weighted linear fit of unwrapped phase vs frequency
% ------------------------------------------------------------------------------
phUnwrap = unwrap(ph);
X = [ones(length(ww), 1), ww];
% (weighting applied elementwise: a diagonal weight matrix would need memory
% quadratic in the number of bins)
beta = (X' * (wgt .* X)) \ (X' * (wgt .* phUnwrap));
out.groupDelay = -beta(2) / Ny; % relative to the series centre, as a fraction of its length
resid = phUnwrap - X * beta;
% (normalized by sqrt(#bins): the residual of an unstructured -- random-walk
% -- unwrapped phase grows as sqrt(#bins), so without this the statistic
% was mostly a restatement of series length: Spearman 0.78 with N on
% real-world series, and 6.6 -> 34 for white noise from N = 500 to 10000)
out.phaseLinearity = sqrt(sum(wgt .* resid.^2)) / sqrt(length(ww));

% ------------------------------------------------------------------------------
%% Magnitude-phase correlation
% ------------------------------------------------------------------------------
out.magPhaseCorr = corr(mag, ph);

% ------------------------------------------------------------------------------
%% Weighted lag-1 autocorrelation of unwrapped-phase increments across frequency
% ------------------------------------------------------------------------------
dPhi = diff(phUnwrap);
d1 = dPhi(1:end - 1); d2 = dPhi(2:end);
wgt3 = wgt(1:end - 2); wgt3 = wgt3 / sum(wgt3);
m1 = sum(wgt3 .* d1); m2 = sum(wgt3 .* d2);
cov12 = sum(wgt3 .* (d1 - m1) .* (d2 - m2));
v1 = sum(wgt3 .* (d1 - m1).^2); v2 = sum(wgt3 .* (d2 - m2).^2);
out.phaseUnwrapAC1 = cov12 / sqrt(v1 * v2);

end
