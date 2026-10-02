function out = NW_VisibilityGraph(y, meth, maxL)
% NW_VisibilityGraph   Visibility graph analysis of a time series.
%
% Constructs a visibility graph of the time series, with one node per sample,
% and returns statistics on the distribution of the number of links per node
% (the degree). In the natural visibility graph ('norm'), two samples are linked
% if the straight line between them passes above every sample in between. In
% the horizontal visibility graph ('horiz'), they are linked if a horizontal
% line between them passes above every sample in between. The outputs
% summarize the degrees (mode, mean, spread, extremes, and heaviness of the
% upper tail), the entropy of their histogram, fits of Gaussian, exponential
% and power-law curves to that histogram and of an extreme-value distribution
% to the degrees, and the autocorrelation of the sequence of degrees taken in
% time order.
%
% ---INPUTS:
% y, the time series (a column vector)
% meth, the method for constructing the graph (default: 'horiz'):
%           (i) 'norm': the natural visibility definition
%           (ii) 'horiz': uses only horizontal lines to link nodes/datums
% maxL, the maximum number of samples to consider (default: 20000). Only the
%       first maxL points of a longer time series are analyzed (a warning is
%       raised), to bound the computation time. Set to 'full' to analyze the
%       entire time series with no cropping (a warning is raised, but no
%       cropping occurs, if the series exceeds 50000 samples, since computation
%       of the natural visibility graph may be slow). Only the degrees are
%       computed, with no adjacency matrix stored, so memory is not a concern:
%       the horizontal graph takes O(N) time and the natural graph typically
%       well under O(N^2) (a smooth random walk of 20000 samples takes about
%       0.5 s; a worst-case series, with no early stop in the sweep, takes O(N^2)).
%
% ---OUTPUTS:
% modek, propmode: the most common degree, and the proportion of nodes that
%       have it
% meank, mediank, stdk: the mean, median and standard deviation of the degrees
% maxk, mink, rangek, iqrk: the maximum, minimum, range and interquartile
%       range of the degrees
% skewnessk: the skewness of the degrees
% maxonmedian: the maximum degree divided by the median degree
% ol90: the mean of the degrees between the 5th and 95th percentiles, divided
%       by the mean of all degrees
% olu90: how far the mean of the top 5% of degrees lies above the overall
%       mean, in standard deviations of the degrees
% dgaussk_r2, dgaussk_adjr2, dgaussk_rmse, dgaussk_resAC1, dgaussk_resAC2,
% dgaussk_resruns: goodness of fit (R^2, adjusted R^2, root-mean-square error),
%       autocorrelation of the residuals at lags 1 and 2, and a runs test
%       p-value, for a single Gaussian fitted to the histogram of degrees
%       (DN_SimpleFit, with as many bins as the range of the degrees); the
%       root-mean-square error is in units of probability density of the degrees
%       divided by their standard deviation, so it does not depend on the number
%       of nodes
% dexpk_r2, dexpk_adjr2, dexpk_rmse, dexpk_resAC1, dexpk_resAC2,
% dexpk_resruns: the same, for a single exponential fitted to the histogram
% dpowerk_r2, dpowerk_adjr2, dpowerk_rmse, dpowerk_resAC1, dpowerk_resAC2,
% dpowerk_resruns: the same, for a power law fitted to the histogram
% gaussnlogL, expnlogL: the mean negative log-likelihood per node of a Gaussian
%       and of an exponential distribution fitted to the degrees
% evparam1, evparam2, evnlogL: the location and scale parameters of an
%       extreme-value distribution fitted to the degrees, and its mean
%       negative log-likelihood per node
% entropy: the entropy of the histogram of degrees (EN_DistributionEntropy,
%       with square-root binning), in nats
% kac1, kac2, kac3: the autocorrelation of the degree sequence, in time order,
%       at lags 1, 2 and 3
% ktau: the lag at which the autocorrelation of the degree sequence first
%       crosses zero (interpolated)
%
% ---REFERENCES:
% Lacasa, Luque, Ballesteros, Luque and Nuno, "From time series to complex
% networks: The visibility graph", P. Natl. Acad. Sci. USA 105(13), 4972
% (2008).
%
% Luque, Lacasa, Ballesteros and Luque, "Horizontal visibility graphs: Exact
% results for random time series", Phys. Rev. E 80(4), 046103 (2009).

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
%% Preliminaries, check inputs
% ------------------------------------------------------------------------------
N = length(y); % time-series length

if size(y, 2) > size(y, 1), y = y'; end % make sure a column vector
if nargin < 2
	% compute the horizontal visibility graph by default
	meth = 'horiz';
end
if nargin < 3 || isempty(maxL)
	maxL = 20000; % crops time series longer than this maximum length
end

if ischar(maxL) && strcmp(maxL, 'full')
	% No cropping -- but flag potentially slow computations for very long series:
	slowThreshold = 50000;
	if N > slowThreshold && strcmp(meth, 'norm')
		warning(sprintf(['Time series (%u samples) exceeds %u with maxL=''full''; ' ...
						 'visibility graph computation may be slow'], N, slowThreshold));
	end
elseif N > maxL % too long: crop
	% ++BF changed on 8/3/2010 to reduce down to first maxL samples. In future,
	% could alter to take different subsets, or set a maximum distance range
	% allowed to make a link (using sparse), etc.
	warning(sprintf(['Time series (%u > %u) is too long for visibility graph...' ...
					 ' Analyzing the first %u samples'], N, maxL, maxL));
	y = y(1:maxL);
	N = length(y); % new time-series length
end

y = y - min(y); % adjust so that minimum of y is at zero

% ------------------------------------------------------------------------------
%% Compute the visibility graph:
% ------------------------------------------------------------------------------
switch meth
	case 'norm'
		% Natural visibility graph degrees by a forward sweep from each node i that keeps the
		% largest slope seen so far: node j is visible from i iff its slope (y(j)-y(i))/(j-i)
		% strictly exceeds that of every node in between. No adjacency matrix is stored.
		% Once the running maximum slope m is positive, nodes beyond distance (max(y)-y(i))/m
		% would have to lie above max(y), so the scan stops early.
		k = zeros(1, N);
		ymax = max(y);
		for i = 1:(N - 1)
			yi = y(i);
			m = -Inf; % largest slope from i seen so far
			jlim = N; % last node that can still be visible
			j = i + 1;
			while j <= jlim
				sj = (y(j) - yi) / (j - i);
				if sj > m
					m = sj;
					k(j) = k(j) + 1;
					k(i) = k(i) + 1;
					if m > 0
						jlim = min(N, i + floor((ymax - yi) / m) + 1);
					end
				end
				j = j + 1;
			end
		end

	case 'horiz'
		% Horizontal visibility graph degrees in O(N) using a monotone stack of
		% non-increasing values. Nodes j < i are linked iff every node between is
		% strictly below min(y(j), y(i)).
		k = zeros(1, N);
		stack = zeros(N, 1);
		sp = 0;
		for i = 1:N
			yi = y(i);
			poppedEqual = false;
			% every stacked node at or below y(i) sees i (the first node ahead of it that is >= it)
			while sp > 0 && y(stack(sp)) <= yi
				j = stack(sp);
				sp = sp - 1;
				k(j) = k(j) + 1;
				k(i) = k(i) + 1;
				poppedEqual = (y(j) == yi);
			end
			% the nearest higher node behind i also sees i, unless an equal-valued node
			% (popped just above) blocks it
			if sp > 0 && ~poppedEqual
				j = stack(sp);
				k(j) = k(j) + 1;
				k(i) = k(i) + 1;
			end
			sp = sp + 1;
			stack(sp) = i;
		end
	otherwise
		error('Unknown visibility graph method ''%s''', meth);
end

% ------------------------------------------------------------------------------
%%% Statistics on the output
% ------------------------------------------------------------------------------

% -------------------------------------------------------------------------------
%% Degree distribution: basic statistics
% -------------------------------------------------------------------------------
out.modek = mode(k); % mode of degree distribution
out.propmode = sum(k == mode(k)) / length(k); % proportion of nodes at the modal degree
out.meank = mean(k); % mean number of links per node
out.mediank = median(k); % median number of links per node
out.stdk = std(k); % std of k
out.maxk = max(k); % maximum degree
out.mink = min(k); % minimum degree
out.rangek = range(k); % range of degree distribution
out.iqrk = iqr(k); % interquartile range of degree distribution
out.skewnessk = skewness(k); % skewness of degree distribution
out.maxonmedian = max(k) / median(k); % max on median (indicator of outlier)
out.ol90 = mean(k(k >= quantile(k, 0.05) & k <= quantile(k, 0.95))) / mean(k);
out.olu90 = (mean(k(k >= quantile(k, 0.95))) - mean(k)) / std(k); % top 5% of points are
% how far from mean (in std units)?

% ------------------------------------------------------------------------------
%% Fit distributions to degree distribution
% ------------------------------------------------------------------------------
% (1) Gauss1: Gaussian fit to degree distribution
try
	dgaussout = DN_SimpleFit(k, 'gauss1', range(k)); % range(k)-bin single gaussian fit
catch emsg
	warning(sprintf('Error fitting gaussian distribution to data:\n%s', emsg.message))
	dgaussout = NaN;
end

if ~isstruct(dgaussout) && isnan(dgaussout)
	out.dgaussk_r2 = NaN;
	out.dgaussk_adjr2 = NaN;
	out.dgaussk_rmse = NaN;
	out.dgaussk_resAC1 = NaN;
	out.dgaussk_resAC2 = NaN;
	out.dgaussk_resruns = NaN;
else
	out.dgaussk_r2 = dgaussout.r2; % rsquared
	out.dgaussk_adjr2 = dgaussout.adjr2; % degrees of freedom-adjusted rsqured
	out.dgaussk_rmse = dgaussout.rmse;  % root mean square error
	out.dgaussk_resAC1 = dgaussout.resAC1; % autocorrelation of residuals at lag 1
	out.dgaussk_resAC2 = dgaussout.resAC2; % autocorrelation of residuals at lag 2
	out.dgaussk_resruns = dgaussout.resruns; % runs test on residuals -- outputs p-value
end

% (2) Exponential1: Exponential fit to degree distribution
try
	dexpout = DN_SimpleFit(k, 'exp1', range(k)); % range(k)-bin single exponential fit
catch emsg
	warning(sprintf('Error fitting exponential distribution to data:\n%s', emsg.message))
	dexpout = NaN;
end
if ~isstruct(dexpout) && isnan(dexpout)
	out.dexpk_r2 = NaN;
	out.dexpk_adjr2 = NaN;
	out.dexpk_rmse = NaN;
	out.dexpk_resAC1 = NaN;
	out.dexpk_resAC2 = NaN;
	out.dexpk_resruns = NaN;
else
	out.dexpk_r2 = dexpout.r2; % rsquared
	out.dexpk_adjr2 = dexpout.adjr2; % degrees of freedom-adjusted rsqured
	out.dexpk_rmse = dexpout.rmse;  % root mean square error
	out.dexpk_resAC1 = dexpout.resAC1; % autocorrelation of residuals at lag 1
	out.dexpk_resAC2 = dexpout.resAC2; % autocorrelation of residuals at lag 2
	out.dexpk_resruns = dexpout.resruns; % runs test on residuals -- outputs p-value
end

% (3) Power1: Power-law fit to degree distribution
try
	dpowerout = DN_SimpleFit(k, 'power1', range(k)); % range(k)-bin single power law fit
catch emsg
	warning(sprintf('Error fitting power-law distribution to data:\n%s', emsg.message))
	dpowerout = NaN;
end
if ~isstruct(dpowerout) && isnan(dpowerout)
	out.dpowerk_r2 = NaN;
	out.dpowerk_adjr2 = NaN;
	out.dpowerk_rmse = NaN;
	out.dpowerk_resAC1 = NaN;
	out.dpowerk_resAC2 = NaN;
	out.dpowerk_resruns = NaN;
else
	out.dpowerk_r2 = dpowerout.r2; % rsquared
	out.dpowerk_adjr2 = dpowerout.adjr2; % degrees of freedom-adjusted rsqured
	out.dpowerk_rmse = dpowerout.rmse;  % root mean square error
	out.dpowerk_resAC1 = dpowerout.resAC1; % autocorrelation of residuals at lag 1
	out.dpowerk_resAC2 = dpowerout.resAC2; % autocorrelation of residuals at lag 2
	out.dpowerk_resruns = dpowerout.resruns; % runs test on residuals -- outputs p-value
end

% ------------------------------------------------------------------------------
%% Using likelihood now:
% ------------------------------------------------------------------------------
% normlike/explike/evlike return the negative log-likelihood *summed* over the
% degree sequence, which is extensive: it grows in direct proportion to the
% number of nodes (and hence to the time-series length), swamping any
% distributional information. Report the mean NLL per node instead, which is
% the intensive quantity these fields were always meant to capture.
numNodes = length(k);

% Gaussian
out.gaussnlogL = normlike([mean(k), std(k)], k) / numNodes;

% Exp
out.expnlogL = explike(mean(k), k) / numNodes;

% Extreme Value Distribution
paramhat = evfit(k);
out.evparam1 = paramhat(1);
out.evparam2 = paramhat(2);
out.evnlogL = evlike(paramhat, k) / numNodes;

% ------------------------------------------------------------------------------
%% Entropy of distribution:
% ------------------------------------------------------------------------------
out.entropy = EN_DistributionEntropy(k, 'hist', 'sqrt');

% Autocorrelations:
out.kac1 = CO_AutoCorr(k, 1, 'Fourier');
out.kac2 = CO_AutoCorr(k, 2, 'Fourier');
out.kac3 = CO_AutoCorr(k, 3, 'Fourier');
out.ktau = CO_FirstCrossing(k, 'ac', 0, 'continuous');

end
