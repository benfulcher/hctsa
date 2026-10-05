function out = CO_CompareMinAMI(y, binMethod, numBins)
% CO_CompareMinAMI   Variability of the first minimum of the automutual information.
%
% Finds the first minimum of the automutual information (AMI), as estimated by
% CO_HistogramAMI, for each number of histogram bins in numBins, up to a maximum
% lag of round(N/2) (set to this maximum if no minimum is found). The function
% returns statistics on the resulting set of first-minimum lags, which measure how
% sensitive the AMI timescale is to the coarse-graining of the time series.
%
% ---INPUTS:
% y, the input time series
% binMethod, the method for estimating mutual information (the meth input to
%            CO_HistogramAMI): 'even', 'std1', 'std2' or 'quantiles'
% numBins, the numbers of bins to compare over (a scalar or a vector; default 10)
%
% ---OUTPUTS:
% min, max, range, median, mean, std, the minimum, maximum, range, median, mean and
%       standard deviation of the first-minimum lags,
% nunique, the number of unique first-minimum lags,
% mode, the most common first-minimum lag,
% modef, the proportion of the bin numbers that give that most common lag,
% conv4, the mean first-minimum lag for the last five bin numbers,
% nprompeaks, the number of prominent peaks (local maxima) of the first-minimum lag as
%       a function of the number of bins: peaks that rise at least 5% of the range of
%       the lags above the higher of the valleys on either side of them, so a small
%       fluctuation does not add a peak.

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
%% Check inputs and set defaults:
% ------------------------------------------------------------------------------
% Default number of bins
if nargin < 3,
	numBins = 10;
end
% ------------------------------------------------------------------------------

doPlot = 0; % plot outputs to figure
N = length(y); % time-series length

% Range of time lags, tau, to consider
%   (although loop usually broken before this maximum)
tauRange = (0:1:round(N / 2));
numTaus = length(tauRange);

% Range of bin numbers to consider
numBinsRange = length(numBins);
amiMins = zeros(numBinsRange, 1);

% Calculate automutual information
for i = 1:numBinsRange % vary over number of bins in histogram
	amis = zeros(numTaus, 1);
	for j = 1:numTaus % vary over time lags, tau
		amis(j) = CO_HistogramAMI(y, tauRange(j), binMethod, numBins(i));
		if (j > 2) && ((amis(j) - amis(j - 1)) * (amis(j - 1) - amis(j - 2)) < 0)
			amiMins(i) = tauRange(j - 1);
			break
		end
	end
	if amiMins(i) == 0
		amiMins(i) = tauRange(end);
	end
end

% Plot:
if doPlot
	figure('color', 'w');
	plot(numBins, amiMins, 'o-k');
end

% -------------------------------------------------------------------------------
% Basic statistics
% -------------------------------------------------------------------------------
out.min = min(amiMins);
out.max = max(amiMins);
out.range = range(amiMins);
out.median = median(amiMins);
out.mean = mean(amiMins);
out.std = std(amiMins);

% Unique values, mode
out.nunique = length(unique(amiMins));
[out.mode, out.modef] = mode(amiMins);
out.modef = out.modef / numBinsRange;

% Converged value?
out.conv4 = mean(amiMins(max(1, end - 4):end));

% -------------------------------------------------------------------------------
% Look for peaks (local maxima)
% -------------------------------------------------------------------------------
% inspired by curious result of periodic maxima for periodic signal with
% bin size... ('quantiles', [2:80])
% Only prominent peaks are counted: a count of every local maximum, or of those
% above a fixed height such as the mean plus one standard deviation, changes with
% each small fluctuation of the curve
out.nprompeaks = SUB_NumProminentPeaks(amiMins, 0.05 * out.range);

end

% -------------------------------------------------------------------------------
function numPeaks = SUB_NumProminentPeaks(x, minProminence)
	% The number of local maxima of x (interior points higher than both neighbors,
	% a flat top counting once) whose prominence is at least minProminence: how far
	% the peak rises above the higher of the lowest values reached on each side
	% before meeting a higher value (or the end of the series).
	x = x(:);
	x = x([true; diff(x) ~= 0]); % merge runs of equal values
	numPeaks = 0;
	for i = 2:length(x) - 1
		if x(i) > x(i - 1) && x(i) > x(i + 1)
			iHigherLeft = find(x(1:i - 1) > x(i), 1, 'last');
			iHigherRight = i + find(x(i + 1:end) > x(i), 1, 'first');
			if isempty(iHigherLeft), iHigherLeft = 1; end
			if isempty(iHigherRight), iHigherRight = length(x); end
			valley = max(min(x(iHigherLeft:i)), min(x(i:iHigherRight)));
			if x(i) - valley >= minProminence
				numPeaks = numPeaks + 1;
			end
		end
	end
end
