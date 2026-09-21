function out = NL_TakensEstimator(y, Nref, rad, past, embedParams, randomSeed)
% NL_TakensEstimator   Taken's estimator for correlation dimension.
%
% cf. "Detecting strange attractors in turbulence", F. Takens.
% Lect. Notes Math. 898 p366 (1981)
%
% Takens' maximum-likelihood estimator of the correlation dimension at an
% upper length scale eup = rad standard deviations of y:
%   D_T = 1 / mean( ln(eup / r_ij) ),
% the mean taken over all pairs (i,j) of delay vectors with max-norm distance
% r_ij < eup, excluding pairs closer in time than the Theiler window, past.
% Computed natively (KD-tree range search at the one radius needed). This
% replaced running TISEAN's d2 (correlation sums over every dimension 1:m
% and every radius, all reference points) followed by c2t, and reading one
% number off the result -- the same estimator, since c2t's
% D_T(r) = C(r) / int_0^r C(r')/r' dr' reduces to the expression above, but
% c2t evaluates it from d2's logarithmically-binned correlation sum at the
% first bin above eup, whereas here it is evaluated exactly at eup from the
% pair distances themselves; values therefore differ slightly from the
% TISEAN-based implementation (which itself was not numerically identical to
% the TSTOOL takens_estimator used before that). Kantz & Schreiber's
% recommendation of half a standard deviation for the length scale is used
% the same way in NL_d2.m's takens05.
%
% ---INPUTS:
% y, the input time series
% Nref, the number of reference points (can be -1 to use all points)
% rad, the upper length scale to read off the dimension estimate, in standard
%       deviations of y (cf. TSTOOL's rad, a proportion of attractor size)
% past, the Theiler window
% embedParams, the embedding parameters for BF_Embed, in the form {tau,m}
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed
%               (relevant if an embedding-dimension method requiring
%               randomization is used)
%
% ---OUTPUT: the Taken's estimator of the correlation dimension, d2 (NaN if
%           no pair of delay vectors lies within the length scale, or all such
%           pairs are exact duplicates, e.g., heavily quantized data).

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

% 2) Upper length scale (standard deviations of y) at which to read off the
% dimension estimate from the correlation-sum data:
if nargin < 3 || isempty(rad)
	rad = 0.05;
end

% 3) Theiler window
if nargin < 4 || isempty(past)
	past = 1; % just exclude current point
end
if (past > 0) && (past < 1)
	past = floor(N * past); % specify a fraction of the time series length...
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

% ------------------------------------------------------------------------------
%% Embed
% ------------------------------------------------------------------------------
Y = BF_Embed(y, embedParams{1}, embedParams{2}, false, randomSeed);
if isscalar(Y) && isnan(Y)
	warning('Could not embed this time series with these embedding parameters');
	out = NaN; return
end
Nemb = size(Y, 1);

% Reference points: the first Nref delay vectors (as TISEAN's d2 -N did), or all:
if Nref == -1 || Nref >= Nemb
	refIdx = (1:Nemb)';
else
	refIdx = (1:Nref)';
end

eup = rad * std(y); % upper length scale, in data units
if ~(eup > 0)
	out = NaN; return % constant series
end

% ------------------------------------------------------------------------------
%% Accumulate sum of ln(eup/r_ij) over pairs with r_ij < eup (max norm),
%% outside the Theiler window, over reference points in chunks (a
%% low-dimensional attractor can have O(N^2) pairs within eup, so never hold
%% them all at once)
% ------------------------------------------------------------------------------
searcher = KDTreeSearcher(Y, 'Distance', 'chebychev');
chunkSize = 500;
sumLog = 0;
numPairs = 0;
for c = 1:chunkSize:length(refIdx)
	theRefs = refIdx(c:min(c + chunkSize - 1, length(refIdx)));
	[idxCell, distCell] = rangesearch(searcher, Y(theRefs, :), eup);
	for k = 1:length(theRefs)
		keep = abs(idxCell{k} - theRefs(k)) > past; % outside the Theiler window (and not self)
		d = distCell{k}(keep);
		d = d(d > 0); % exact duplicates carry no length-scale information (ln -> Inf)
		sumLog = sumLog + sum(log(eup ./ d));
		numPairs = numPairs + length(d);
	end
end

if numPairs == 0
	warning('No pairs within %g standard deviations of each other to estimate a correlation dimension from', rad);
	out = NaN; return
end
% (No minimum pair count beyond that. Note that for high embedding
% dimensions of noise-like series no pair may fall within eup at all, giving
% NaN here (~1/3 of Empirical1000 series for the m = 8 and m = 10 variants);
% the TISEAN-based implementation returned a constant 14.3 in that situation
% -- 1/ln of d2's radius-bin ratio, i.e. every pair in a single bin -- which
% was an artifact, not an estimate. Where both are defined they agree to
% ~2% (median), Spearman 0.95-0.99 across Empirical1000.)

out = numPairs / sumLog; % Takens' estimator: 1 / mean(ln(eup/r))

end
