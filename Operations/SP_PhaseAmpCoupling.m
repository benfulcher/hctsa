function out = SP_PhaseAmpCoupling(y, nBands, maxN, nPhaseBins)
% SP_PhaseAmpCoupling   Cross-frequency phase-amplitude coupling
%
% Splits the (one-sided, DC- and Nyquist-excluded) spectrum into nBands
% equal-width frequency bands -- the same equal-band convention
% SP_Summaries uses for its band-power fields -- and, for every pair of
% bands i < j, asks whether the instantaneous phase of the slower band i
% modulates the instantaneous amplitude envelope of the faster band j.
% This is the standard "phase-amplitude coupling" (PAC) measure from the
% neuroscience literature (e.g. theta phase modulating gamma amplitude in
% EEG), generalized here to arbitrary equal-width bands rather than
% fixed, domain-specific ones.
%
% Each band's instantaneous phase/amplitude is obtained directly from an
% FFT-domain analytic signal: zeroing every bin outside the band and
% doubling the surviving positive-frequency bins before an inverse FFT
% gives the band-limited analytic signal in one step (the same
% zero-negative-frequencies construction behind the textbook FFT-based
% Hilbert transform), with no Signal Processing Toolbox dependency.
%
% Coupling within each band pair is quantified by Tort et al.'s
% modulation index (MI): the faster band's amplitude envelope is binned
% by the slower band's instantaneous phase, and MI is the
% Kullback-Leibler divergence of that phase-binned mean-amplitude
% distribution from uniform, normalized to [0,1] by its maximum
% (log(nPhaseBins)) -- 0 for amplitude independent of phase, higher as
% amplitude becomes concentrated at a preferred phase.
%
% cf. SP_Bicoherence, which detects quadratic *phase* coupling between
% frequency triplets (f1, f2, f1+f2); this instead detects *amplitude*
% modulation of one band by the phase of another, a distinct
% cross-frequency mechanism that a bicoherence/bispectrum analysis does
% not target.
%
% ---INPUTS:
% y, the input time series
%
% nBands, the number of equal-width frequency bands to split the
%         spectrum into (default: 5, matching SP_Summaries' 5-band
%         split). Phase-amplitude pairs are formed from every pair of
%         bands i < j (phase from the slower band, amplitude from the
%         faster), giving nchoosek(nBands,2) pairs.
%
% maxN, the maximum number of samples to consider (cf. SP_Bicoherence's
%       maxN); longer series are cropped to their first maxN points.
%       Can be 'full' to disable cropping (default).
%
% nPhaseBins, the number of phase bins used to estimate each band pair's
%             modulation index (default: 18, i.e. 20-degree bins, the
%             standard choice from Tort et al. 2010).
%
% ---OUTPUTS:
% maxMI, the maximum modulation index across all band pairs -- the
%       comodulogram peak, i.e. whether *any* band pair shows real
%       coupling (cf. SP_Bicoherence's maxBic). Chosen over the mean or
%       standard deviation of MI across pairs after checking all three
%       on Empirical1000: meanMI and stdMI were both near-duplicates of
%       maxMI (r=0.96 and r=0.99 respectively) *and* mechanically
%       diluted by nBands -- doubling nBands from 5 to 10 roughly halved
%       meanMI (more near-zero pairs enter the average) while maxMI
%       barely moved (x0.97), making it the only one of the three whose
%       meaning doesn't depend on this parameter choice.
%
% entropyMI, the normalized Shannon entropy of the MI values across
%       pairs (0 = coupling concentrated in a single band pair, 1 =
%       uniformly diffuse across all pairs) -- cf. SP_Bicoherence's
%       analogous entropy field for the bicoherence surface. The one
%       field found to carry genuinely separate information from maxMI
%       (r=0.08-0.23 on Empirical1000, vs. r>=0.95 for meanMI/stdMI).

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

if nargin < 2 || isempty(nBands)
    nBands = 5;
end
if nargin < 3 || isempty(maxN)
    maxN = 'full';
end
if nargin < 4 || isempty(nPhaseBins)
    nPhaseBins = 18;
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
%% Equal-width frequency bands (DC and Nyquist bins excluded, as in
%% SP_SummariesPhase -- neither carries meaningful oscillatory phase)
% ------------------------------------------------------------------------------
halfN = floor(N / 2) + 1; % one-sided bins: DC (1) up to Nyquist (halfN, if N even)
isNyquistBin = (mod(N, 2) == 0);
if isNyquistBin
    usableBins = 2:(halfN - 1);
else
    usableBins = 2:halfN;
end

minBinsPerBand = 2; % need >=2 FFT bins per band for a non-degenerate analytic signal
if length(usableBins) < nBands * minBinsPerBand
    warning('Time series too short (N = %u) for %u usable frequency bands', N, nBands);
    out = NaN; return
end

edges = round(linspace(1, length(usableBins) + 1, nBands + 1));
bandBins = cell(nBands, 1);
for b = 1:nBands
    bandBins{b} = usableBins(edges(b):edges(b + 1) - 1);
end

% ------------------------------------------------------------------------------
%% Per-band analytic signal (FFT-domain Hilbert trick, band-limited in one step)
% ------------------------------------------------------------------------------
Y = fft(y);
phaseBand = cell(nBands, 1);
ampBand = cell(nBands, 1);
for b = 1:nBands
    Yb = zeros(N, 1);
    Yb(bandBins{b}) = 2 * Y(bandBins{b});
    zb = ifft(Yb);
    phaseBand{b} = angle(zb);
    ampBand{b} = abs(zb);
end

% ------------------------------------------------------------------------------
%% Modulation index (Tort et al. 2010) for every phase(i)-amplitude(j) pair, i < j
% ------------------------------------------------------------------------------
phaseEdges = linspace(-pi, pi, nPhaseBins + 1);
Hmax = log(nPhaseBins);
numPairs = nchoosek(nBands, 2);
MI = zeros(numPairs, 1);
k = 0;
for i = 1:nBands - 1
    for j = i + 1:nBands
        k = k + 1;
        phi = phaseBand{i};
        A = ampBand{j};
        [~, binIdx] = histc(phi, phaseEdges);
        binIdx(binIdx == nPhaseBins + 1) = nPhaseBins; % phi == pi edge case
        meanAmp = accumarray(binIdx, A, [nPhaseBins, 1], @mean, 0);
        p = meanAmp / sum(meanAmp);
        p = p(p > 0);
        H = -sum(p .* log(p));
        MI(k) = (Hmax - H) / Hmax;
    end
end

if all(~isfinite(MI))
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Summary statistics across band pairs
% ------------------------------------------------------------------------------
out.maxMI = max(MI);

% Normalized Shannon entropy of the MI values across pairs (0 = coupling
% concentrated in one pair, 1 = uniformly diffuse):
p = MI(MI > 0);
if isempty(p)
    out.entropyMI = NaN;
else
    p = p / sum(p);
    out.entropyMI = -sum(p .* log(p)) / log(numPairs);
end

end
