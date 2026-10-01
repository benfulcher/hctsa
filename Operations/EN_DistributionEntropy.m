function out = EN_DistributionEntropy(y, histOrKS, numBins, olremp)
% EN_DistributionEntropy   Entropy of the distribution of values in a time series.
%
% Estimates the entropy of the distribution of a data vector, ignoring the order
% of values in time. The distribution is estimated either with a histogram or as
% a kernel-smoothed density (using ksdensity from Matlab's Statistics Toolbox,
% evaluated on a grid of 200 points spanning the 0.1%-99.9% quantile range plus a
% 10% margin). The entropy is -sum(p.*log(p./w)), where p is the probability in
% each cell and w the cell width, so it estimates a differential entropy (in
% nats) and depends on the scale of the data. For the histogram estimate, a
% Miller-Madow correction for the finite sample, (number of nonempty bins - 1)/
% (2*length(y)), is added.
%
% An optional additional parameter can be used to remove a proportion of the most
% extreme values at both ends of the distribution as an initial preprocessing.
%
% ---INPUTS:
% y, the input time series
% histOrKS, 'hist' for a histogram, or 'ks' for a kernel-smoothed density
%    (default: 'hist')
% numBins, for 'hist': either a positive integer, giving the number of
%        equal-width bins, or the name of a rule for choosing the bin width
%        passed to histcounts ('auto', 'fd', 'sqrt', 'sturges', ...);
%        for 'ks': a positive real number, the width parameter for ksdensity, or
%        empty for the default (automatically chosen) width, which is optimal for
%        a Gaussian distribution
%        (default: 10)
% olremp [optional], the proportion of values to remove at both extremes (by
%        quantile; e.g., olremp = 0.01 keeps only the middle 98% of the data; 0
%        keeps all data). This parameter ought to be less than 0.5, which keeps
%        none of the data. If olremp is nonzero, the output is the difference in
%        entropy from removing the outliers (full data minus trimmed data).
%        (default: 0)
%
% ---OUTPUTS:
% a scalar: the entropy estimate (in nats), or, if olremp is nonzero, the
% entropy of the full time series minus that of the trimmed time series.
% NaN if everything is removed by the trimming, or if the 'ks' grid range is
% degenerate (near-constant series).

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

doPlot = 0; % plot outputs to figure

% ------------------------------------------------------------------------------
%% Check inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(histOrKS)
	histOrKS = 'hist'; % use histogram by default
end
if nargin < 3 % (can be empty for default width for ksdensity)
	numBins = 10; % use 10 bins
end
if nargin < 4
	olremp = 0;
end

% ------------------------------------------------------------------------------
% (1) Remove outliers?
% ------------------------------------------------------------------------------
if olremp ~= 0
	yHat = y(y >= quantile(y, olremp) & y <= quantile(y, 1 - olremp));
	if isempty(yHat)
		% removed the entire time series?!
		% shouldn't be possible for good values of olremp with equality
		% in the above inequalities
		out = NaN; return
	else
		% Return the difference in entropy from removing outliers
		out = EN_DistributionEntropy(y, histOrKS, numBins) - ...
				EN_DistributionEntropy(yHat, histOrKS, numBins);
		return
	end
end

% ------------------------------------------------------------------------------
% (2) Form the histogram
% ------------------------------------------------------------------------------
switch histOrKS
	case 'hist' % Use histogram to calculate pdf
		if isnumeric(numBins)
			[px, binEdges] = histcounts(y, numBins, 'Normalization', 'probability');
		else
			[px, binEdges] = histcounts(y, 'BinMethod', numBins, 'Normalization', 'probability');
		end
		% Compute bin centers:
		xr = mean([binEdges(1:end - 1); binEdges(2:end)]);
		% Compute bin widths:
		binWidths = diff(binEdges);

	case 'ks' % Use ksdensity to calculate pdf
		% Evaluate on an explicit, length-stable grid.
		%
		% ksdensity's default grid spans the observed data range, and the range
		% of a sample grows with N (as ~sqrt(2*log(N)) for Gaussian data), so
		% the grid -- and with it the log(binWidth) term in the entropy sum
		% below -- widened with time-series length regardless of the underlying
		% distribution. Anchoring the grid to extreme *quantiles* instead fixes
		% that: quantiles are consistent estimators, so the interval converges
		% as N grows rather than expanding.
		numGridPts = 200;
		lo = quantile(y, 0.001);
		hi = quantile(y, 0.999);
		if ~(hi > lo) % degenerate (near-constant) input
			out = NaN; return
		end
		pad = 0.1 * (hi - lo); % a little headroom beyond the quantile range
		xGrid = linspace(lo - pad, hi + pad, numGridPts);
		if isempty(numBins)
			[px, xr] = ksdensity(y, xGrid, 'function', 'pdf'); % selects optimal width
		else
			% NB: a *fixed* absolute bandwidth makes the density estimate
			% inconsistent -- for consistency the bandwidth must shrink with
			% sample size (Silverman: h ~ N^(-1/5)) -- so the smoothness of the
			% estimated density, and hence its entropy, drifts with N whatever
			% the underlying distribution. Measured eta^2 against length was
			% 0.56-0.97 for fixed bandwidths against 0.08 for the automatic
			% selection, so the fixed-bandwidth variants are no longer
			% registered as hctsa features. The option is kept for callers who
			% want a specific smoothing scale.
			[px, xr] = ksdensity(y, xGrid, 'width', numBins, 'function', 'pdf'); % uses specified width
		end
		binWidths = ones(1, length(px)) * (xr(2) - xr(1));
		% ksdensity returns a *density* evaluated on a grid, whereas the entropy
		% sum below (shared with the 'hist' branch) expects probability mass per
		% cell. Using the raw density there left sum(px) ~= 1, so the result was
		% neither a discrete nor a differential entropy. Convert to probability
		% mass and renormalize (the grid truncates a little tail mass).
		px = px .* binWidths;
		px = px / sum(px);

	otherwise
		error('Unknown distribution method -- specify ''ks'' or ''hist''') % error; must specify 'ks' or 'hist'
end

if doPlot
	figure('color', 'w'); box('on');
	plot(xr, px, 'k')
end

% ------------------------------------------------------------------------------
% (3) Compute the entropy sum and return it as output
% ------------------------------------------------------------------------------
% 0*log0 = 0:
% -sum(p.*log(p./binWidth)) = H_discrete + mean log binWidth, i.e. the standard
% discretized differential entropy.
out = -sum(px(px > 0) .* log(px(px > 0) ./ binWidths(px > 0)));

% Miller-Madow correction to the discrete part, for the histogram estimator
% only: the plug-in entropy is biased low by (M-1)/(2n) for M occupied bins and
% n samples, which made the finer binnings (20 and 50 bins) track time-series
% length more than the distribution. Kernel density estimates are smoothed
% rather than plug-in, so the discrete-bin correction does not apply to them.
if strcmp(histOrKS, 'hist')
	out = out + (sum(px > 0) - 1) / (2 * length(y));
end

end
