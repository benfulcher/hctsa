function out = SP_SpectralTimeFreq(y, numWindows)
% SP_SpectralTimeFreq  Time-varying spectral statistics from a spectrogram.
%
% SP_Summaries computes statistics from a single, static spectral estimate of
% the whole time series. This function instead divides the series into
% overlapping windows and tracks how the spectral content changes across
% them, using Matlab's Signal Processing Toolbox:
%
% (i) Spectral kurtosis: for each frequency bin, the kurtosis of that bin's
%     power across all windows. High values flag a frequency band whose
%     energy is concentrated in occasional bursts rather than spread evenly
%     over time (e.g., a transient, impulsive fault).
%
% (ii) Instantaneous spectral entropy: the Shannon entropy of the power
%      spectrum computed separately in each window, giving one entropy value
%      per window. Variation in this sequence flags a time series whose
%      spectral character is not stationary.
%
% ---INPUTS:
%
% y, the input time series
%
% numWindows, the target number of overlapping windows (50% overlap) to
%             divide the series into (default: 20)

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
	error('Time series (N=%u) too short for a %u-window spectrogram', Ny, numWindows);
end

% ------------------------------------------------------------------------------
% Spectral kurtosis: kurtosis across windows, per frequency bin
% ------------------------------------------------------------------------------
% MATLAB version compatibility: spectralKurtosis output signatures vary across
% releases (e.g., older releases can return fewer outputs).
try
	[kurt, spread, centroid, thresh, fout] = spectralKurtosis(y, Fs, ...
				'Window', window, 'OverlapLength', noverlap, ...
				'Scaled', false, 'ConfidenceLevel', 0.95);
catch emsg
	if contains(emsg.message, 'Too many output arguments')
		% Older MATLAB signatures: keep core outputs and mark unavailable ones NaN.
		[kurt, fout] = spectralKurtosis(y, Fs, 'Window', window, 'OverlapLength', noverlap);
        spread = NaN;
        centroid = NaN;
        thresh = NaN;
	else
		rethrow(emsg)
	end
end

out.sk_max = max(kurt);
out.sk_mean = mean(kurt);
out.sk_std = std(kurt);
out.sk_range = max(kurt) - min(kurt);
out.sk_fracAboveThresh = mean(kurt > thresh); % fraction of frequencies with non-Gaussian, bursty behavior
[~, i_max] = max(kurt);
out.sk_freqAtMax = 2 * pi * fout(i_max); % angular frequency, matching SP_Summaries convention
if all(isnan(spread))
	out.sk_meanSpread = NaN;
else
	out.sk_meanSpread = mean(spread);
end
if all(isnan(centroid))
	out.sk_meanCentroid = NaN;
else
	out.sk_meanCentroid = 2 * pi * mean(centroid);
end

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
