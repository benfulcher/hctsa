function out = ST_PeakIntervals(y, minProm)
% ST_PeakIntervals   Prominence and regularity of the peaks of a time series.
%
% Finds the local maxima of the z-scored series whose topographic prominence is
% at least minProm standard deviations (MATLAB's findpeaks with
% 'MinPeakProminence'; a peak's prominence is its height above the higher of the
% two lowest points separating it from any higher peak, or from the series end).
% Peaks are located in time, so this is a time-domain counterpart to the spectral
% peak summaries of SP_Summaries, and it ignores small wiggles that
% ST_LocalExtrema or a zero-crossing count would register. It then summarizes the
% typical prominence of the peaks and how regular the spacing between successive
% peaks is: the inter-peak intervals of a periodic or strongly quasi-periodic
% series are nearly equal (CV near 0), those of a point process are broadly
% distributed, and a series with bursts of peaks has serially correlated intervals.
%
% Intervals are in samples. Samples at the two ends of the series cannot be
% peaks (findpeaks does not report them).
%
% Requires the Signal Processing Toolbox (findpeaks).
%
% ---INPUTS:
% y, the input time series (z-scored internally, so minProm is in standard
%    deviations)
% minProm, the minimum peak prominence, in standard deviations (default: 1)
%
% ---OUTPUTS:
% meanProm, the mean prominence of the detected peaks (NaN if there are none)
% cvInt, the coefficient of variation of the inter-peak intervals, std/mean (NaN with fewer than 3 intervals)
% acInt1, the lag-1 correlation between successive inter-peak intervals (NaN with fewer than 5
%         intervals, or if the intervals are all equal)
% All fields are NaN for constant, non-finite, or very short (N < 20) series.
%
% ---NOTES:
% For a z-scored pure sinusoid the mean prominence is 2*sqrt(2) = 2.83 and
% cvInt is 0. The peak rate (peaks per sample) was computed during development
% but not kept: it is almost perfectly rank-correlated (r > 0.97) with existing
% features (zero-crossing and local-difference measures) at minProm = 0.5, 1 and 2.

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

% findpeaks is in the Signal Processing Toolbox
BF_CheckToolbox('signal_toolbox');

y = y(:);
if nargin < 2 || isempty(minProm)
    minProm = 1;
end

out.meanProm = NaN;
out.cvInt = NaN;
out.acInt1 = NaN;

N = length(y);
if N < 20 || any(~isfinite(y)) || std(y) == 0
    return
end
y = (y - mean(y)) / std(y); % prominence is in units of the series' standard deviation

[~, locs, ~, proms] = findpeaks(y, 'MinPeakProminence', minProm);
nPk = length(locs);

if nPk >= 1
    out.meanProm = mean(proms);
end

ipi = diff(locs); % inter-peak intervals (samples)
if length(ipi) >= 3
    out.cvInt = std(ipi) / mean(ipi);
end
if length(ipi) >= 5
    r = corrcoef(ipi(1:end-1), ipi(2:end));
    out.acInt1 = r(1, 2); % NaN if all intervals are equal
end

end
