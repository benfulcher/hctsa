function out = SP_PhaseFluctuationScaling(y, halfWidthFrac, numWindows, maxN)
% SP_PhaseFluctuationScaling    Multi-scale fluctuation analysis of instantaneous phase
%
% Isolates the dominant oscillatory component of y (the frequency band
% carrying the most spectral power, excluding DC and Nyquist), takes its
% instantaneous phase via the analytic signal, and asks how the
% fluctuation of that phase about its mean rotation rate grows with
% window size N: mean(|dphi(t+N) - dphi(t)|) vs. N, in log-log space (the
% same detrended-fluctuation-analysis logic SC_FluctAnal/SY_ apply to raw
% values, applied here to instantaneous phase instead).
%
% This targets a specific finding of Gupta, Prasad, Singh & Ramaswamy
% (below): for strange nonchaotic dynamics, this curve rises at short
% windows (reflecting local instability -- SNAs and chaotic dynamics look
% identical over short windows) but *flattens* at longer windows, because
% an SNA's largest Lyapunov exponent is negative and its phase dynamics
% are, over long timescales, globally stable. For chaotic dynamics, the
% curve keeps rising at long windows too, because nearby trajectories
% (and hence the phase) never stop diverging. This is a real signature to
% distinguish from a similar-looking but different thing: an ordinary,
% cleanly periodic oscillation also flattens quickly, but its short-window
% slope is close to zero as well (there's no local instability to begin
% with) -- it's the combination of an appreciable short-window rise *and*
% a later flattening that is closer to the paper's SNA signature, not
% flattening alone.
%
% The original method uses empirical mode decomposition (EMD) to isolate
% the dominant intrinsic mode. This implementation instead bandpasses
% around the single dominant non-edge FFT peak (the same band-limited
% analytic-signal trick SP_PhaseAmpCoupling uses), which is simpler,
% parameter-transparent, and avoids EMD's known mode-mixing/boundary
% sensitivities -- at the cost of not being a literal reimplementation.
%
% cf. K. Gupta, A. Prasad, H.P. Singh, R. Ramaswamy, "Analytical signal
% analysis of strange nonchaotic dynamics", Phys. Rev. E 77, 046220 (2008).
% DOI: 10.1103/PhysRevE.77.046220
%
% ---INPUTS:
%
% y, the input time series
%
% halfWidthFrac, half-width of the frequency band around the dominant
%                peak, as a fraction of the usable (DC- and
%                Nyquist-excluded) one-sided spectrum (default: 0.01)
%
% numWindows, the number of log-spaced window sizes N to evaluate, split
%             evenly into a short-window half and a long-window half for
%             two separate linear fits in log-log space (default: 16)
%
% maxN, the maximum number of samples to consider (cf. NL_RQA; default:
%       10000)
%
% ---OUTPUTS: the short-window and long-window log-log slopes
% (slope_short, slope_long), their difference (slope_diff = slope_short -
% slope_long, positive and appreciable for the SNA-like flattening
% signature above), and the mean rotation frequency of the isolated
% dominant component (meanFreq, in cycles per sample -- correlates r=0.90
% with the existing SP_Summaries_*_maxw fields on the Bonn EEG dataset,
% which capture a related "location of the dominant spectral peak"
% concept via a plain periodogram rather than this operation's isolated
% band; kept anyway since it's a free byproduct of the computation above
% and the two aren't identical).

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
%% Check inputs
% ------------------------------------------------------------------------------
N = length(y);
y = y(:);

if nargin < 2 || isempty(halfWidthFrac)
    halfWidthFrac = 0.01;
end
if nargin < 3 || isempty(numWindows)
    numWindows = 16;
end
if nargin < 4 || isempty(maxN)
    maxN = 10000;
end

if N > maxN
    warning('Time series (%u samples) exceeds maxN = %u; analyzing the first %u samples', N, maxN, maxN);
    y = y(1:maxN);
    N = maxN;
end

% ------------------------------------------------------------------------------
%% Locate the dominant non-edge spectral peak (DC and Nyquist bins
%% excluded, as in SP_SummariesPhase/SP_PhaseAmpCoupling)
% ------------------------------------------------------------------------------
halfN = floor(N / 2) + 1;
isNyquistBin = (mod(N, 2) == 0);
if isNyquistBin
    usableBins = 2:(halfN - 1);
else
    usableBins = 2:halfN;
end

halfWidthBins = max(2, round(halfWidthFrac * length(usableBins)));
if length(usableBins) < 4 * halfWidthBins
    warning('Time series too short (N = %u) for a well-defined dominant frequency band', N);
    out = NaN; return
end

Y = fft(y);
power = abs(Y(usableBins)).^2;
[~, peakIdxRel] = max(power);
peakBin = usableBins(peakIdxRel);

bandLo = max(usableBins(1), peakBin - halfWidthBins);
bandHi = min(usableBins(end), peakBin + halfWidthBins);
bandBins = bandLo:bandHi;

% ------------------------------------------------------------------------------
%% Band-limited analytic signal (FFT-domain Hilbert trick) and its
%% instantaneous phase
% ------------------------------------------------------------------------------
Yband = zeros(N, 1);
Yband(bandBins) = 2 * Y(bandBins);
z = ifft(Yband);
phi = unwrap(angle(z));

if any(~isfinite(phi))
    warning('Non-finite instantaneous phase (degenerate band-limited signal?)');
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Detrend: remove the mean rotation rate to leave the phase fluctuation
% ------------------------------------------------------------------------------
tIdx = (1:N)';
A = [tIdx, ones(N, 1)];
coef = A \ phi;
out.meanFreq = coef(1) / (2 * pi); % cycles per sample
dphi = phi - A * coef;

% ------------------------------------------------------------------------------
%% Multi-scale fluctuation analysis of dphi: mean(|dphi(t+w) - dphi(t)|)
%% vs. window size w, in log-log space, split into a short-window and a
%% long-window linear fit
% ------------------------------------------------------------------------------
wMin = 2;
wMax = floor(N / 4);
if wMax <= wMin * 4
    warning('Time series too short (N = %u) for a meaningful multi-scale phase-fluctuation analysis', N);
    out.slope_short = NaN;
    out.slope_long = NaN;
    out.slope_diff = NaN;
    return
end

windows = unique(round(logspace(log10(wMin), log10(wMax), numWindows)));
if length(windows) < 6
    out.slope_short = NaN;
    out.slope_long = NaN;
    out.slope_diff = NaN;
    return
end

meanAbsDiff = zeros(length(windows), 1);
for k = 1:length(windows)
    w = windows(k);
    d = abs(dphi(1 + w:end) - dphi(1:end - w));
    meanAbsDiff(k) = mean(d);
end

logW = log10(windows(:));
logD = log10(meanAbsDiff);
valid = isfinite(logD) & (meanAbsDiff > 0);
if sum(valid) < 6
    out.slope_short = NaN;
    out.slope_long = NaN;
    out.slope_diff = NaN;
    return
end
logW = logW(valid);
logD = logD(valid);

splitPoint = ceil(length(logW) / 2);
shortIdx = 1:splitPoint;
longIdx = (splitPoint):length(logW); % one point of overlap anchors the two fits together

out.slope_short = SUB_fitSlope(logW(shortIdx), logD(shortIdx));
out.slope_long = SUB_fitSlope(logW(longIdx), logD(longIdx));
out.slope_diff = out.slope_short - out.slope_long;

% ------------------------------------------------------------------------------
function slope = SUB_fitSlope(x, yv)
    if length(x) < 2
        slope = NaN; return
    end
    Afit = [x(:), ones(length(x), 1)];
    coefFit = Afit \ yv(:);
    slope = coefFit(1);
end
% ------------------------------------------------------------------------------

end
