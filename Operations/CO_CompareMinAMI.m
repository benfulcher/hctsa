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
% nlocmax, the number of local maxima of the first-minimum lag, as a function of
%       the number of bins, that lie more than one standard deviation above the mean.

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
% local maxima above 1*std from mean
% inspired by curious result of periodic maxima for periodic signal with
% bin size... ('quantiles', [2:80])
loc_extr = intersect(find(diff(amiMins(1:end - 1)) > 0), BF_SignChange(diff(amiMins(1:end - 1)), 1)) + 1;
big_loc_extr = intersect(find(amiMins > out.mean + out.std), loc_extr);
out.nlocmax = length(big_loc_extr);

end
