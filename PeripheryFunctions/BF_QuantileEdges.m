function edges = BF_QuantileEdges(y, numBins)
% BF_QuantileEdges   Histogram bin edges at the quantiles of the data (equiprobable bins).
%
% Returns the edges, to pass to histcounts, of numBins bins with about the same
% number of values in each: the quantiles of y at 0, 1/numBins, ..., 1 (as computed
% by quantile). Repeated quantiles (tied values) are merged, giving fewer bins. The
% interior edges are lowered, and the end edges widened, by a tiny fraction (1e-9) of
% the range of the data, since a quantile often equals a data value, which a
% rounding error in the last bit could otherwise move between neighboring bins;
% values on an edge always fall just above it, in the upper bin.
%
% ---INPUTS:
% y, the data vector (NaNs are ignored)
% numBins, the number of bins
%
% ---OUTPUTS:
% edges, a row vector of increasing bin edges.

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

y = y(~isnan(y));
edges = unique(quantile(y, linspace(0, 1, numBins + 1)));
if length(edges) == 1 % constant data: one bin of unit width
	edges = edges + [-0.5, 0.5];
	return
end
tol = 1e-9 * (edges(end) - edges(1));
edges(1:end - 1) = edges(1:end - 1) - tol; % lower all but the last edge,
edges(end) = edges(end) + tol; % and widen the last

end
