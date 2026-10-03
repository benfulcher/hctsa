function out = NL_FNN(y, tau, maxm, theilerWin, justBest, bestp, escapeFactor)
% NL_FNN   How the fraction of false nearest neighbors falls as the embedding dimension of the series increases.
%
% Uses the false_nearest routine from the TISEAN package for nonlinear
% time-series analysis. For each embedding dimension m = 1,...,maxm (delay tau),
% every point's nearest neighbor (in the maximum norm) is found, excluding points
% within the Theiler window in time. The neighbor is false if, when one more
% coordinate is added, the two points move apart by more than escapeFactor times
% their original distance. For a deterministic system the fraction of false
% neighbors falls toward 0 once the embedding dimension is large enough; for
% noise it stays high. Neighbors farther apart than the standard deviation of
% the data divided by the escape factor are skipped.
%
% The TISEAN routines are run in the command line using 'system' commands in
% MATLAB, and require that TISEAN is installed and compiled, and able to be
% executed in the command line.
%
% ---INPUTS:
% y, the input time series
% tau, the time delay (a number of samples, or 'ac' for the first zero-crossing
%      of the autocorrelation function, 'ac1e' for the floor of its first 1/e
%      crossing, or 'mi' for the smaller of the first minimum of the Kraskov
%      automutual information and the 'ac1e' delay; see BF_GetTau; default: 1)
% maxm, the maximum embedding dimension
% theilerWin, the Theiler window: {'ac', k} for k times the first zero-crossing
%             of the autocorrelation function, or a number of samples (see
%             BF_TheilerWindow; default: {'ac',1})
% justBest, if 1 just outputs a scalar estimate of embedding dimension: the first
%           dimension at which the fraction of false nearest neighbors is below
%           bestp (default: 1)
% bestp, only used if justBest==1 -- the fnn threshold for picking an embedding
%        dimension (default: 0.4)
% escapeFactor [opt], the neighbor-distance escape factor (TISEAN's '-f',
%                its R_tol-like false-neighbor threshold: a candidate neighbor
%                is "false" if its distance grows by more than this factor
%                after one more embedding dimension). Default is TISEAN's own
%                internal default of 2.0 (matches false_nearest's behavior when
%                -f is omitted). Note this is notably stricter than the
%                analogous 'th' parameter of MS_fnn/NL_MS_fnn, whose default
%                is 5 -- at matched escapeFactor values the two methods give
%                closely comparable per-dimension false-neighbor profiles.
%
% ---OUTPUTS: (if justBest is 1, a scalar embedding dimension; otherwise a
% structure with the fields below, where i = 1,...,maxm is the embedding dimension)
% pfnn_<i>: the fraction of false nearest neighbors in an i-dimensional embedding
% nHood2_<i>: the typical (root-mean-square) distance from a point to its nearest
%             neighbor in an i-dimensional embedding
% minpfnn, meanpfnn, stdpfnn: minimum, mean and standard deviation of the
%             fraction of false nearest neighbors across dimensions
% maxnHood2, meannHood2: maximum and mean of nHood2 across dimensions
% firstunder09, firstunder08, firstunder07, firstunder06, firstunder05,
%             firstunder04, firstunder03, firstunder02, firstunder01,
%             firstunder005: the first embedding dimension at which the fraction
%             of false nearest neighbors is below 90%, 80%, ..., 10%, 5% (maxm + 1
%             if it never is)
% max1stepchange: the largest absolute change in the fraction between
%             consecutive embedding dimensions
% mdrop: the mean change in the fraction per added dimension
% pdrop: minus the mean sign of the change, i.e. the fraction of decreasing steps
%        minus the fraction of increasing steps
%
% ---REFERENCES:
% R. Hegger, H. Kantz and T. Schreiber, "Practical implementation of nonlinear
% time series methods: The TISEAN package", Chaos 9(2), 413 (1999).
%
% ---NOTES:
% TISEAN is available at
% http://www.mpipks-dresden.mpg.de/~tisean/Tisean_3.0.1/index.html and the
% false_nearest documentation at
% http://www.mpipks-dresden.mpg.de/~tisean/TISEAN_2.1/docs/docs_c/false_nearest.html
%
% The fourth column of TISEAN's output (nHood2) is the square root of the mean
% squared nearest-neighbor distance.
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

doPlot = false; % can turn on to see plotted summaries

% ------------------------------------------------------------------------------
%% Check inputs / set defaults
% ------------------------------------------------------------------------------
if nargin < 1
	error('Input a time series')
end
N = length(y);
if N < 10
	warning('Time series (N=%u) too short for fnn', N);
	out = NaN; return
end

if nargin < 2 || isempty(tau)
	tau = 1; % time delay
end
if ischar(tau) && strcmp(tau, 'ac1e')
	% Adaptive delay: see BF_GetTau
	tau = BF_GetTau(y, tau);
	if isnan(tau)
		out = NaN; return
	end
end
if strcmp(tau, 'ac')
	tau = CO_FirstCrossing(y, 'ac', 0, 'discrete'); % first zero-crossing of autocorrelation function
elseif strcmp(tau, 'mi')
	tau = BF_GetTau(y, 'mi'); % min(first Kraskov AMI minimum, 1/e ACF time)
end
if isnan(tau)
	out = NaN; return
end

% Maximum embedding dimension:
if nargin < 3
	maxm = 10;
end

% Theiler window:
if nargin < 4 || isempty(theilerWin)
	theilerWin = {'ac', 1};
end
theilerWin = BF_TheilerWindow(y, theilerWin, N);
if isnan(theilerWin) % the autocorrelation function never crosses zero
	warning('No autocorrelation zero-crossing to set the Theiler window')
	out = NaN; return
end

% Just return best dimension:
if nargin < 5 || isempty(justBest)
	justBest = true; % just return the best embedding dimension
end

% How to return the best embedding dimension:
if nargin < 6
	bestp = 0.4; % stop when under 40% false nearest neighbors
end

% Escape factor (TISEAN's '-f'; omitted here to preserve false_nearest's own
% internal default of 2.0 unless explicitly overridden):
if nargin < 7 || isempty(escapeFactor)
	escapeFactor = [];
end

% ------------------------------------------------------------------------------
%% Write the file
% ------------------------------------------------------------------------------
filePath = BF_WriteTempFile(y);
if TISEANVerbose()
	fprintf(1, 'Wrote the input time series (N = %u) to the temporary file ''%s'' for TISEAN.\n', length(y), filePath);
end

% ------------------------------------------------------------------------------
%% Run the TISEAN code, false_nearest
% ------------------------------------------------------------------------------
if isempty(escapeFactor)
	tisean_command = sprintf('false_nearest -d%u -m1 -M1,%u -t%u -V0 %s', tau, maxm, theilerWin, filePath);
else
	tisean_command = sprintf('false_nearest -d%u -m1 -M1,%u -t%u -f%g -V0 %s', tau, maxm, theilerWin, escapeFactor, filePath);
end
% BF_TiseanSystem already error()s if false_nearest couldn't be run at all
% (not installed, or hung and hit its own timeout guard -- verified this
% binary can block indefinitely on a missing/unreadable input file rather
% than erroring); any status it returns here means the binary ran to
% completion, whether cleanly or with its own terse refusal (see below).
[~, res] = BF_TiseanSystem(tisean_command);

% first column: the embedding dimension
% second column: the fraction of false nearest neighbors
% third column: the average size of the neighborhood
% fourth column: the average of the squared size of the neighborhood

% Read TISEAN output:
% (data-dependent failure, not a code error: e.g. a large ACF-based tau on a
% strongly-autocorrelated/near-unit-root series -- such as a random walk --
% can require (maxm-1)*tau points, more than N provides, so TISEAN refuses
% and prints an explanatory message to stdout instead of data rows. Return
% NaN rather than error(), matching the rest of the codebase's convention;
% BF_Embed and its callers already guard for a NaN embedding dimension.)
if isempty(res)
	warning('No output from TISEAN routine false_nearest on the data');
	out = NaN; return
end
data = textscan(res, '%u%f%f%f');

mDim = double(data{1}); % embedding dimension
pNN = data{2}; % fraction of false nearest neighbors
% nHoodSize = data{3}; % average size of the neighbourhood
nHoodSize2 = data{4}; % average squared size of the neighbourhood

% Check that some data exists:
if isempty(mDim) || isempty(pNN) || isempty(nHoodSize2)
	warning('TISEAN false_nearest produced no usable output for this data (e.g. requested embedding dimension/delay too large for the series length)');
	out = NaN; return
end

if doPlot
	f = figure('color', 'w'); box('on'); hold on
	plot(pNN, 'o-k'); plot(nHoodSize2, 'o-r')
	legend('pNN', 'mean squared size of neighbourhood')
	xlabel('Embedding dimension');
end

% ------------------------------------------------------------------------------
% Output(s)
% ------------------------------------------------------------------------------
if justBest
	% We just want a scalar to choose the embedding with
	out = firstunderf(bestp, mDim, pNN);
	return
end

% Output all of them
for i = 1:maxm
	if i <= length(mDim)
		out.(sprintf('pfnn_%u', i)) = pNN(i); % proportion of false nearest neighbors
		out.(sprintf('nHood2_%u', i)) = nHoodSize2(i); % mean squared size of neighbourhood
	else
		% Not enough points found to estimate at this dimension
		out.(sprintf('pfnn_%u', i)) = NaN;
		out.(sprintf('nHood2_%u', i)) = NaN;
	end
end

% pNN summaries:
out.minpfnn = min(pNN); % minimum
out.meanpfnn = mean(pNN); % mean
out.stdpfnn = std(pNN); % standard deviation

% nHood2 summaries:
out.maxnHood2 = max(nHoodSize2); % maximum
out.meannHood2 = mean(nHoodSize2); % mean

% Find embedding dimension for the first time p goes under x%
out.firstunder09 = firstunderf(0.9, mDim, pNN);   % 90%
out.firstunder08 = firstunderf(0.8, mDim, pNN);   % 80%
out.firstunder07 = firstunderf(0.7, mDim, pNN);   % 70%
out.firstunder06 = firstunderf(0.6, mDim, pNN);   % 60%
out.firstunder05 = firstunderf(0.5, mDim, pNN);   % 50%
out.firstunder04 = firstunderf(0.4, mDim, pNN);   % 40%
out.firstunder03 = firstunderf(0.3, mDim, pNN);   % 30%
out.firstunder02 = firstunderf(0.2, mDim, pNN);   % 20%
out.firstunder01 = firstunderf(0.1, mDim, pNN);   % 10%
out.firstunder005 = firstunderf(0.05, mDim, pNN); % 5%

% Maximum step-wise change across p
out.max1stepchange = max(abs(diff(pNN)));

% Curve shape across dimensions (mean signed per-step change, and the
% proportion of m -> m+1 steps for which pfnn decreased):
out.mdrop = mean(diff(pNN));
out.pdrop = -mean(sign(diff(pNN)));

% ------------------------------------------------------------------------------
function firsti = firstunderf(x, m, p)
	%% Find m for the first time p goes under x%
	firsti = m(find(p < x, 1, 'first'));
	if isempty(firsti)
		firsti = m(length(m)) + 1;
	end
end

end
