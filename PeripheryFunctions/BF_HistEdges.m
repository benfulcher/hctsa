function edges = BF_HistEdges(y, binRule, limits)
% BF_HistEdges   Equal-width histogram bin edges from an explicit bin-count rule.
%
% Returns the bin edges, to pass to histcounts, for a histogram of equal-width
% bins spanning the data, with the number of bins either given or set by a rule
% written out as a formula (histcounts' own BinMethod rules round the bin width to
% 'nice' values, so the bins depend on the implementation).
%
% The bins span [min(y), max(y)] (or limits) with width (max - min)/numBins. The interior edges
% are lowered, and the end edges widened, by a tiny fraction (1e-6) of a bin width,
% so that values lying exactly on an edge of the ideal grid (e.g., lattice-valued
% data or counts), which a rounding error in the last bit could otherwise move
% between neighboring bins, always fall just above it, in the upper bin.
%
% ---INPUTS:
% y, the data vector (NaNs are ignored)
% binRule, the number of bins (a positive integer), or a rule for it, for n values:
%       'sqrt': ceil(sqrt(n)),
%       'sturges': ceil(log2(n) + 1),
%       'fd': Freedman-Diaconis, ceil(range/(2*IQR*n^(-1/3))) (Sturges if IQR = 0),
%       'auto': the larger of the Sturges and Freedman-Diaconis numbers of bins.
%       (default: 'auto'). A rule's number of bins is at most n, and at least 1.
% limits, [lower, upper], the interval to span, instead of the range of the data
%       (optional; the number of bins from a rule still depends on the data)
%
% ---OUTPUTS:
% edges, a row vector of numBins + 1 increasing bin edges.

% ------------------------------------------------------------------------------
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

if nargin < 2 || isempty(binRule)
	binRule = 'auto';
end
y = y(~isnan(y));
n = length(y);
dataRange = max(y) - min(y); % (for the rules)
if nargin >= 3 && ~isempty(limits)
	lo = limits(1);
	hi = limits(2);
else
	lo = min(y);
	hi = max(y);
end

% ------------------------------------------------------------------------------
% Number of bins
% ------------------------------------------------------------------------------
if isnumeric(binRule)
	numBins = binRule;
else
	numSturges = ceil(log2(n) + 1);
	if dataRange > 0
		fdWidth = 2 * iqr(y) * n^(-1/3);
	else
		fdWidth = 0;
	end
	if fdWidth > 0
		numFD = ceil(dataRange / fdWidth);
	else
		numFD = numSturges;
	end
	switch binRule
		case 'sqrt'
			numBins = ceil(sqrt(n));
		case 'sturges'
			numBins = numSturges;
		case 'fd'
			numBins = numFD;
		case 'auto'
			numBins = max(numSturges, numFD);
		otherwise
			error('Unknown bin rule ''%s''', binRule)
	end
	numBins = max(1, min(numBins, n));
end

% ------------------------------------------------------------------------------
% Bin edges
% ------------------------------------------------------------------------------
if hi == lo % constant data: one bin of unit width
	edges = lo + [-0.5, 0.5];
	return
end
binWidth = (hi - lo) / numBins;
tol = 1e-6 * binWidth;
edges = lo + (0:numBins) * binWidth - tol; % interior edges, lowered by tol
edges(1) = lo - tol;
edges(end) = hi + tol;

end
