function out = DN_HistogramMode(y,numBins)
% DN_HistogramMode      Mode of a data vector.
%
% Measures the mode of the data vector using histograms with a given number
% of bins.
%
%---INPUTS:
%
% y, the input data vector
%
% numBins, the number of bins to use in the histogram.

% ------------------------------------------------------------------------------
% Copyright (C) 2017, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

if nargin < 2
    numBins = 'auto'; % ceil(sqrt(length(y)));
end

% Compute the histogram from the data:
if isnumeric(numBins)
    [N,binEdges] = histcounts(y,numBins);
elseif ischar(numBins)
    [N,binEdges] = histcounts(y,'BinMethod',numBins);
else
    error('Unknown format for numBins');
end

% Compute bin centers from bin edges:
binCenters = mean([binEdges(1:end-1); binEdges(2:end)]);

% Mean position of maximums (if multiple):
out = mean(binCenters(N == max(N)));

end
