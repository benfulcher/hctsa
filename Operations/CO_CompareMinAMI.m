function out = CO_CompareMinAMI(y,binMethod,numBins)
% CO_CompareMinAMI	Variability in first minimum of automutual information
%
% Finds the first minimum of the automutual information by various different
% estimation methods, and sees how this varies over different coarse-grainings
% of the time series.
%
% The function returns a set of statistics on the set of first minimums of the
% automutual information function obtained over a range of the number of bins
% used in the histogram estimation, when specifying 'numBins' as a vector
%
%---INPUTS:
% y, the input time series
%
% binMethod, the method for estimating mutual information (input to CO_HistogramAMI)
%
% numBins, the number of bins for the AMI estimation to compare over (can be a
%           scalar or vector)
%
% Outputs include the minimum, maximum, range, number of unique values, and the
% position and periodicity of peaks in the set of automutual information
% minimums.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
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
%	(although loop usually broken before this maximum)
tauRange = (0:1:round(N/2));
numTaus = length(tauRange);

% Range of bin numbers to consider
numBinsRange = length(numBins);
amiMins = zeros(numBinsRange,1);

% Calculate automutual information
for i = 1:numBinsRange % vary over number of bins in histogram
    amis = zeros(numTaus,1);
    for j = 1:numTaus % vary over time lags, tau
        amis(j) = CO_HistogramAMI(y,tauRange(j),binMethod,numBins(i));
        if (j > 2) && ((amis(j)-amis(j-1))*(amis(j-1)-amis(j-2)) < 0)
            amiMins(i) = tauRange(j-1);
            break
        end
    end
    if amiMins(i) == 0
		amiMins(i) = tauRange(end);
	end
end

% Plot:
if doPlot
    figure('color','w');
    plot(numBins,amiMins,'o-k');
end

%-------------------------------------------------------------------------------
% Basic statistics
%-------------------------------------------------------------------------------
out.min = min(amiMins);
out.max = max(amiMins);
out.range = range(amiMins);
out.median = median(amiMins);
out.mean = mean(amiMins);
out.std = std(amiMins);

% Unique values, mode
out.nunique = length(unique(amiMins));
[out.mode, out.modef] = mode(amiMins);
out.modef = out.modef/numBinsRange;

% Converged value?
out.conv4 = mean(amiMins(end-4:end));

%-------------------------------------------------------------------------------
% Look for peaks (local maxima)
%-------------------------------------------------------------------------------
% local maxima above 1*std from mean
% inspired by curious result of periodic maxima for periodic signal with
% bin size... ('quantiles', [2:80])
loc_extr = intersect(find(diff(amiMins(1:end-1)) > 0), BF_sgnchange(diff(amiMins(1:end-1)),1)) + 1;
big_loc_extr = intersect(find(amiMins > out.mean+out.std),loc_extr);
out.nlocmax = length(big_loc_extr);

end
