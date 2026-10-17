function out = NL_TakensEstimator(y, Nref, c0, past, embedParams, randomSeed)
% NL_TakensEstimator   Takens' estimator of the correlation dimension, at the length scale enclosing a fixed fraction of pairs.
%
% Estimates the correlation dimension, D, of a time-delay embedding of y with Takens'
% estimator, at an upper length scale eup set by the data: eup is the c0 quantile of the
% max-norm distances between pairs of delay vectors (pairs closer in time than the Theiler
% window, and repeated points, are excluded). With r_ij the n pair distances below eup,
%   D = (n - 1) / sum( ln(eup / r_ij) ).
% Under the estimator's model, C(r) proportional to r^D below eup, the ln(eup / r_ij) are
% independent exponential variables with rate D given eup (Takens 1985; with eup an order
% statistic of the distances, Hill 1975), so their sum is Gamma(n, D) and (n - 1)/sum is the
% unbiased estimate of D (the maximum-likelihood form n/sum is biased upward by n/(n - 1);
% Borovkova et al. 1999).
%
% Setting eup by a fraction of pairs, rather than as a fixed radius, keeps it at a small
% scale that always holds the same share of pairs, whatever the embedding dimension and the
% spread of the embedded points, and its target (the scale below which a fraction c0 of
% pairs lie) does not depend on the length of the series. (At a fixed radius, noise-like
% series embedded in several dimensions have no pairs at all at small scales.) Takens'
% estimator is a reasonable dimension estimate when evaluated at the upper end of the
% scaling range (Hegger, Kantz & Schreiber 1999), which for a low-dimensional attractor lies
% at small scales.
%
% Details:
% - Pairs of delay vectors closer than tol = 1e-12*(max(y) - min(y)) are repeated points:
%   they carry no length-scale information, and a positive fraction of them makes the
%   estimator inconsistent (Borovkova et al. 1999), so they are left out.
% - Distances closer than tol are equal. If the values of y lie on a lattice, i.e., all
%   differences between distinct values are within 1% of integer multiples of a common step
%   q (the smallest difference divided by the smallest k = 1, ..., 100 that fits, well above
%   tol, with q fitted to the differences by least squares; quantized data, also when stored
%   with rounding), distances closer than q/2 are equal, since finer differences are not
%   resolved by the data. Equal distances are merged into levels after sorting (not snapped
%   to a grid), each level taking its smallest value; pairs at the level of eup lie at eup.
%   The tolerances are relative to the data, so shifting, rescaling or reflecting y leaves
%   the output unchanged, apart from distances within round-off of tol.
% - The output is NaN if fewer than 3 distinct distances lie below eup: the data then do
%   not resolve the length scale (eup lies within two resolution steps of the smallest
%   distance), and the estimate is set by the lattice of the values rather than by the
%   dynamics.
% - The output is NaN if the estimate exceeds 2*m for an m-dimensional embedding: a set with
%   a scaling region below eup has a dimension of at most m, and finite-sample and
%   finite-scale effects keep estimates for sets with a scaling region well below 2*m.
%   Larger values arise when the correlation sum rises much faster than r^m below eup, as
%   when the pairs outside the Theiler window are all far apart (a Theiler window that is
%   a large fraction of a nonstationary series), so there is no scaling region at that
%   scale.
% - The output is NaN if c0 times the number of pairs that are not repeated points is
%   below 33, the smallest number of pairs for which the 95% sampling interval of the
%   estimate (for independent pairs) spans less than a factor of two.
% - The pairs are (i,j) with i among the first Nref delay vectors and j > i + past.
%
% ---INPUTS:
% y, the input time series
% Nref, the number of reference points (the first Nref delay vectors; -1 uses all points;
%       default: -1)
% c0, the fraction of pairs (other than repeated points) that lie below the upper length
%     scale, eup
%     (default: 0.01)
% past, the Theiler window: {'ac', k} for k times the first zero-crossing of the
%       autocorrelation function, or a number of samples (see BF_TheilerWindow; default:
%       {'ac', 1})
% embedParams, the embedding parameters for BF_Embed, in the form {tau,m}, where tau
%              is a number of samples or 'ac', 'ac1e' or 'mi' (see BF_GetTau) and m
%              is a number of dimensions, or 'fnn' or {'fnn',th} (false nearest
%              neighbors at threshold th, at most floor(2*log10(N)); see BF_Embed)
%              (default: {'ac','fnn'})
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed (relevant if
%             an embedding-dimension method requiring randomization is used)
%
% ---OUTPUTS:
% A scalar: the estimate of the correlation dimension (NaN if the series is constant or
% cannot be embedded, the Theiler window cannot be set, c0 times the number of pairs that
% are not repeated points is below 33, fewer than 3 distinct distances lie below eup, or
% the estimate exceeds twice the embedding dimension; see Details).
%
% ---REFERENCES:
% F. Takens, "On the numerical determination of the dimension of an attractor", in
% B.L.J. Braaksma, H.W. Broer and F. Takens (eds.), Dynamical Systems and
% Bifurcations (Groningen, 1984), Lecture Notes in Mathematics 1125, 99-106,
% Springer, Berlin (1985). DOI: 10.1007/BFb0075637
% B.M. Hill, "A simple general approach to inference about the tail of a distribution",
% Ann. Statist. 3(5), 1163-1174 (1975). DOI: 10.1214/aos/1176343247
% S. Borovkova, R. Burton and H. Dehling, "Consistency of the Takens estimator for the
% correlation dimension", Ann. Appl. Probab. 9(2), 376-390 (1999).
% DOI: 10.1214/aoap/1029962747
% R. Hegger, H. Kantz and T. Schreiber, "Practical implementation of nonlinear time series
% methods: The TISEAN package", Chaos 9(2), 413-435 (1999). DOI: 10.1063/1.166424

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
%% Check inputs
% ------------------------------------------------------------------------------
N = length(y); % time-series length

% 1) Nref
if nargin < 2 || isempty(Nref)
	Nref = -1; % use all points
end

% 2) Fraction of pairs below the upper length scale
if nargin < 3 || isempty(c0)
	c0 = 0.01;
end
if ~(isscalar(c0) && c0 > 0 && c0 < 1)
	error('c0, the fraction of pairs below the upper length scale, must be in (0,1)');
end

% 3) Theiler window
if nargin < 4 || isempty(past)
	past = {'ac', 1};
end
past = BF_TheilerWindow(y, past, N);
if isnan(past) % the autocorrelation function never crosses zero
	warning('No autocorrelation zero-crossing to set the Theiler window')
	out = NaN; return
end

% 4) Embedding parameters
if nargin < 5 || isempty(embedParams)
	embedParams = {'ac', 'fnn'};
	fprintf(1, 'Using default time-delay embedding using autocorrelation and fnn-mar\n');
else
	if length(embedParams) ~= 2
		error('Embedding parameters are incorrectly formatted, we need {tau,m}')
	end
end

% 5) randomSeed: how to treat the randomization
if nargin < 6
	randomSeed = []; % default
end

minPairs = 33; % fewest pairs below eup for an estimate

% ------------------------------------------------------------------------------
%% Embed
% ------------------------------------------------------------------------------
Y = BF_Embed(y, embedParams{1}, embedParams{2}, false, randomSeed);
if isscalar(Y) && isnan(Y)
	warning('Could not embed this time series with these embedding parameters');
	out = NaN; return
end
Nemb = size(Y, 1);
if Nref == -1 || Nref >= Nemb
	nRef = Nemb;
else
	nRef = Nref;
end

% Pairs closer than tol = 1e-12 of the range of y are repeated points, left out
span = max(y) - min(y); % no max-norm distance exceeds the range of y
if ~(span > 0)
	out = NaN; return % constant series
end
tol = 1e-12 * span;
% Quantized values: if the differences between distinct values are all within 1% of
% integer multiples of a common step q (the smallest difference divided by the smallest
% k = 1, ..., 100 that fits, well above tol), distances closer than q/2 are not resolved
% by the data
gaps = diff(unique(y(:)));
gaps = gaps(gaps > tol);
if isempty(gaps)
	out = NaN; return % constant to within tol
end
q = 0;
for k = 1:100
	step = min(gaps) / k;
	if step <= 1000*tol, break; end
	mult = floor(gaps/step + 0.5);
	if all(abs(gaps - step*mult) <= 0.01*step)
		q = sum(gaps.*mult) / sum(mult.^2);
		break
	end
end
mergeTol = max(tol, q/2); % distances closer than this are equal

% ------------------------------------------------------------------------------
%% Pass 1: histogram of the distances, to locate the c0 quantile
% ------------------------------------------------------------------------------
numBins = 200;
edges = tol * logspace(0, log10(2*span/tol + 1), numBins + 1);
counts = zeros(1, numBins);
for c = 1:SUB_chunkSize(Nemb):nRef
	d = SUB_distances(Y, c, min(c + SUB_chunkSize(Nemb) - 1, nRef), past, tol, Inf);
	counts = counts + histcounts(d, edges);
end
numPos = sum(counts); % pairs that are not repeated points
K = ceil(c0 * numPos); % the c0 quantile is the K-th smallest distance
if K < minPairs
	warning('Too few pairs (%u below the %g quantile) to estimate a dimension', K, c0);
	out = NaN; return
end

% ------------------------------------------------------------------------------
%% Pass 2: collect the distances below the bin edge beyond the quantile
% ------------------------------------------------------------------------------
b = find(cumsum(counts) >= K, 1, 'first');
rMax = edges(b + 1);
p = [];
for c = 1:SUB_chunkSize(Nemb):nRef
	p = [p; SUB_distances(Y, c, min(c + SUB_chunkSize(Nemb) - 1, nRef), past, tol, rMax)]; %#ok<AGROW>
end
p = sort(p);

% Merge distances closer than mergeTol into levels (each takes its smallest value); eup
% is the level of the K-th distance
level = [0; cumsum(diff(p) > mergeTol)];
value = p([true; diff(level) > 0]);
k = level(K); % the number of levels below eup
if k < 3
	warning('Fewer than 3 distinct pair distances below the upper length scale');
	out = NaN; return % the data are too discrete at this scale
end

% ------------------------------------------------------------------------------
%% The estimate
% ------------------------------------------------------------------------------
eup = value(k + 1); % (levels are numbered from 0)
r = value(level(level < k) + 1); % the pairs below eup, at the values of their levels
out = (numel(r) - 1) / sum(log(eup ./ r));
m = size(Y, 2);
if out > 2*m
	warning('Dimension estimate %g exceeds twice the embedding dimension %u', out, m);
	out = NaN; % no scaling region at this scale
end

end
% ------------------------------------------------------------------------------
function s = SUB_chunkSize(Nemb)
	% Reference points per block (about 2.5 million distances at a time)
	s = max(1, floor(2.5e6 / Nemb));
end
% ------------------------------------------------------------------------------
function d = SUB_distances(Y, i1, i2, past, tol, rMax)
	% Max-norm distances above tol and below rMax between reference points i1:i2 and
	% the delay vectors j > i + past
	Nemb = size(Y, 1);
	j1 = i1 + past + 1;
	if j1 > Nemb
		d = zeros(0, 1); return
	end
	ii = (i1:i2)';
	jj = j1:Nemb;
	D = abs(Y(ii, 1) - Y(jj, 1)');
	for k = 2:size(Y, 2)
		D = max(D, abs(Y(ii, k) - Y(jj, k)'));
	end
	keep = (jj > ii + past) & (D > tol) & (D < rMax);
	d = D(keep);
end
