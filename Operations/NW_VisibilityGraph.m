function out = NW_VisibilityGraph(y, meth, maxL)
% NW_VisibilityGraph    Visibility graph analysis of a time series.
%
% Constructs a visibility graph of the time series and returns various
% statistics on the properties of the resulting network.
%
% cf.: "From time series to complex networks: The visibility graph"
% Lacasa, Lucas and Luque, Bartolo and Ballesteros, Fernando and Luque, Jordi
% and Nuno, Juan Carlos P. Natl. Acad. Sci. USA. 105(13) 4972 (2008)
%
% "Horizontal visibility graphs: Exact results for random time series"
% Luque, B. and Lacasa, L. and Ballesteros, F. and Luque, J.
% Phys. Rev. E. 80(4) 046103 (2009)
%
% ---INPUTS:
%
% y, the time series (a column vector)
%
% meth, the method for constructing:
%           (i) 'norm': the normal visibility definition
%           (ii) 'horiz': uses only horizonatal lines to link nodes/datums
%
% maxL, the maximum number of samples to consider. Due to memory constraints,
%               only the first maxL (5000 by default) points of time series are
%               analyzed. Longer time series are reduced to their first maxL
%               samples. Set to 'full' to analyze the entire time series with
%               no cropping (a warning is raised, but no cropping occurs, if
%               the series exceeds 10000 samples, since computation may be slow).
%
% ---OUTPUTS:
%
% Statistics on the degree distribution, including the mode, mean, spread,
% histogram entropy, and fits to gaussian, exponential, and power-law distributions.

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
	maxL = 5000; % crops time series longer than this maximum length
end

if ischar(maxL) && strcmp(maxL, 'full')
	% No cropping -- but flag potentially slow computations for very long series:
	slowThreshold = 10000;
	if N > slowThreshold
		warning(sprintf(['Time series (%u samples) exceeds %u with maxL=''full''; ' ...
						 'visibility graph computation may be slow'], N, slowThreshold));
	end
elseif N > maxL % too long to store in memory
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
		% Normal (natural) visibility graph (Lacasa et al., 2008).
		% Original code by Enyu Zhuang (Zoey), 23/9/13; substantially
		% modified by Ben Fulcher, 26-8-2015; inlined here and reduced
		% from an N x N gradient matrix to a single reused row vector
		% (only the current row i was ever read), 4/8/2026.
		A = zeros(N, N);
		s = zeros(1, N); % gradient from i to each subsequent point j (reused per i)
		for i = 1:(N - 1)
			for j = i + 1:N
				s(j) = (y(j) - y(i)) / (j - i);
				if j == i + 1
					A(i, j) = 1; % always visible to the next point
				else
					% Visible iff the gradient to every intermediate point
					% is lower than the gradient straight to j:
					for k = i + 1:j - 1
						if s(k) >= s(j)
							break;
						elseif k == j - 1
							A(i, j) = 1;
						end
					end
				end
			end
		end
		A = symmetrize(A);

	case 'horiz'
		% Horizontal visibility graph

		% The graph has only O(N) edges, so accumulate (row,col) edge indices
		% and build a sparse adjacency matrix at the end, instead of a dense
		% N x N matrix (an O(N^2) allocation/initialization that's mostly
		% zeros -- e.g., ~200MB at the default maxL=5000 cap):
		yr = flipud(y); % reversed order time series
		rows = zeros(2 * N, 1);
		cols = zeros(2 * N, 1);
		numEdges = 0;

		for i = 1:N
			% Look forward to first blocker, then stop
			if i < N
				nAhead = find(y(i + 1:end) > y(i), 1, 'first');
				if ~isempty(nAhead)
					numEdges = numEdges + 1;
					rows(numEdges) = i;
					cols(numEdges) = i + nAhead;
				end
			end

			% Look back to the first hit, then stop
			if i > 1
				nBack = find(yr(N - i + 2:end) > yr(N - i + 1), 1, 'first');
				if ~isempty(nBack)
					numEdges = numEdges + 1;
					rows(numEdges) = i - nBack;
					cols(numEdges) = i;
				end
			end
		end
		% Every edge added above has row < col (no self-loops, strictly upper
		% triangular), so collapsing any duplicate (row,col) pairs added from
		% both directions back to 0/1 (sparse() sums duplicates) matches the
		% original dense code's idempotent A(i,j)=1 assignment exactly:
		A = double(sparse(rows(1:numEdges), cols(1:numEdges), 1, N, N) > 0);

		% Symmetrize A (safe because A is strictly upper triangular: A and A'
		% have disjoint nonzero support, so no double-counting):
		A = A + A';
	otherwise
		error('Unknown visibility graph method ''%s''', meth);
end

% ------------------------------------------------------------------------------
%%% Statistics on the output
% ------------------------------------------------------------------------------

% -------------------------------------------------------------------------------
%% Degree distribution: basic statistics
% -------------------------------------------------------------------------------
k = sum(A); % the degree distribution
k = full(k);

out.modek = mode(k); % mode of degree distribution
out.propmode = sum(k == mode(k)) / sum(k);
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

% -------------------------------------------------------------------------------
function A = symmetrize(A)
	% Symmetrize an upper triangular matrix:
	At = A';
	lowerT = logical(tril(ones(size(A))));
	A(lowerT) = At(lowerT);
end
