function out = SP_SpectralTimeFreq(y, numWindows)
% SP_SpectralTimeFreq   Time-varying spectral statistics from a spectrogram.
%
% SP_Summaries computes statistics from a single, static spectral estimate of the
% whole time series. This function instead divides the series into overlapping
% windows and tracks how the spectral content changes across them, using Matlab's
% Signal Processing Toolbox:
%
% (i) Spectral kurtosis: for each frequency bin, the kurtosis of that bin's power
%     across all windows. High values flag a frequency band whose energy is
%     concentrated in occasional bursts rather than spread evenly over time
%     (e.g., a transient, impulsive fault).
%
% (ii) Instantaneous spectral entropy: the Shannon entropy of the power spectrum
%      computed separately in each window, giving one entropy value per window.
%      Variation in this sequence flags a time series whose spectral character is
%      not stationary.
%
% Each window is a Hamming window of max(8, round(N/numWindows)) samples with 50%
% overlap, so about 2*numWindows - 1 windows result.
%
% ---INPUTS:
% y, the input time series
% numWindows, sets the window length to N/numWindows samples (at least 8), with
%             50% overlap (default: 20). If fewer than 4 windows fit, all
%             outputs are NaN.
%
% ---OUTPUTS:
% sk_max, sk_mean, sk_std, sk_range: maximum, mean, standard deviation and range,
%         over frequencies, of the spectral kurtosis
% sk_fracAboveThresh, the fraction of frequencies whose spectral kurtosis exceeds
%         the 95% Gaussian-null threshold (non-Gaussian, bursty behavior)
% sk_freqAtMax, the angular frequency (2*pi*f, matching SP_Summaries) at which the
%         spectral kurtosis is largest
% sk_meanSpread, the mean over frequencies of the standard deviation across windows
%         of the power in each frequency bin (the spread output of
%         spectralKurtosis); the power is |FFT|^2/(0.5*sum(window)^2), so it
%         scales with the variance of the series and as 1/(window length)
% sk_meanCentroid, 2*pi times the mean over frequencies of the mean across windows
%         of the same power (the centroid output of spectralKurtosis, which for
%         unscaled spectral kurtosis is a mean power, not a frequency); not
%         registered as a feature
% se_mean, se_std, se_max, se_min, se_range: mean, standard deviation, maximum,
%         minimum and range, over windows, of the spectral entropy
%
% ---NOTES:
% All outputs are computed directly from the spectrogram, so they do not depend on
% the MATLAB release (on releases without the five-output form of spectralKurtosis,
% sk_meanSpread and sk_meanCentroid used to be NaN). The values equal those of the
% toolbox function's outputs to rounding error.

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
%% Check that a Signal Processing Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('signal_toolbox');

% ------------------------------------------------------------------------------
% Check inputs, set defaults:
% ------------------------------------------------------------------------------
if size(y, 2) > size(y, 1)
	y = y'; % Time series must be a column vector
end
if nargin < 2 || isempty(numWindows)
	numWindows = 20;
end

Ny = length(y); % time-series length
Fs = 1; % sampling frequency (hctsa convention: dimensionless, unit-spaced series)

% ------------------------------------------------------------------------------
% Set the window: sized as a fraction of the series length (rather than a
% fixed number of samples) so behavior scales sensibly across hctsa's very
% different input lengths. Matlab's own default window for these functions
% (rectwin(round(Fs*0.03))) assumes a physical sampling rate in Hz and
% degenerates to a zero-length window at hctsa's Fs = 1 convention.
% ------------------------------------------------------------------------------
winLength = max(8, round(Ny / numWindows));
noverlap = round(winLength / 2);
window = hamming(winLength);

% Need enough windows for the across-window statistics below to be meaningful:
hopLength = winLength - noverlap;
numFrames = floor((Ny - winLength) / hopLength) + 1;
if numFrames < 4
	% Too short for across-window statistics to mean anything: NaN for all outputs
	out = NaN; return
end

% ------------------------------------------------------------------------------
% Spectral kurtosis: kurtosis across windows, per frequency bin
% ------------------------------------------------------------------------------
% Everything is computed directly from the spectrogram, so that the outputs do
% not depend on the MATLAB release (the outputs of spectralKurtosis differ
% between releases, and older ones lack the per-frequency form). The values are
% those of [kurt, spread, centroid, thresh, fout] = spectralKurtosis(y, Fs,
% 'Window', window, 'OverlapLength', noverlap, 'Scaled', false, 'ConfidenceLevel',
% 0.95), which are all functions of P, the power in each frequency bin and
% window, normalized as |FFT|^2/(0.5*sum(window)^2) and halved at zero frequency
% (and at the Nyquist frequency for an even window length), with K windows:
%   kurt = ((K+1)/(K-1)) <P^2>/<P>^2 - 2, means over windows (Antoni 2006);
%   spread = standard deviation of P across windows; and
%   centroid = <P>, the mean of P across windows (named a centroid by MATLAB, but
%              with Scaled = false it is a mean power, not a frequency).
[Sxx, fout] = spectrogram(y, window, noverlap, winLength, Fs);
P = abs(Sxx).^2 / (0.5 * sum(window)^2);
P(1, :) = 0.5 * P(1, :); % zero frequency
if rem(winLength, 2) == 0
	P(end, :) = 0.5 * P(end, :); % Nyquist frequency
end
K = size(P, 2); % number of windows
kurt = ((K + 1) / (K - 1)) * mean(P.^2, 2) ./ mean(P, 2).^2 - 2;
thresh = 2 * sqrt(2) * erfcinv(1 - 0.95) / sqrt(K); % 95% Gaussian-null threshold
spread = std(P, 0, 2);
centroid = mean(P, 2);

out.sk_max = max(kurt);
out.sk_mean = mean(kurt);
out.sk_std = std(kurt);
out.sk_range = max(kurt) - min(kurt);
out.sk_fracAboveThresh = mean(kurt > thresh); % fraction of frequencies with non-Gaussian, bursty behavior
[~, i_max] = max(kurt);
out.sk_freqAtMax = 2 * pi * fout(i_max); % angular frequency, matching SP_Summaries convention
out.sk_meanSpread = mean(spread);
out.sk_meanCentroid = 2 * pi * mean(centroid);

% ------------------------------------------------------------------------------
% Instantaneous spectral entropy: entropy per window, across windows
% ------------------------------------------------------------------------------
se = spectralEntropy(y, Fs, 'Window', window, 'OverlapLength', noverlap);

out.se_mean = mean(se);
out.se_std = std(se);
out.se_max = max(se);
out.se_min = min(se);
out.se_range = max(se) - min(se);

end
