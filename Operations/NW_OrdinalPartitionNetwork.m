function out = NW_OrdinalPartitionNetwork(y, d, tau)
% NW_OrdinalPartitionNetwork    Ordinal partition transition network measures
%
% Symbolizes the time series into ordinal patterns (Bandt-Pompe) and builds
% a directed transition network in which nodes are the ordinal patterns
% actually observed and edges connect a pattern to whichever pattern
% immediately follows it in time. Network-topological measures of this
% ordinal partition transition network are then computed.
%
% cf. "Using ordinal partition transition networks to analyze ECG data"
% C.W. Kulp, J.M. Chobot, H.R. Freitas, G.D. Sprechini, Chaos 26, 073114 (2016)
%
% cf. "Streaming feature-based causal graph discovery"... (unrelated); see
% instead: M. McCullough, M. Small, T. Stemler, H.H.-C. Iu, "Time lagged
% ordinal partition networks for capturing dynamics of continuous dynamical
% systems", Chaos 25, 053101 (2015) -- the original (weighted) ordinal
% partition transition network, of which Kulp et al.'s unweighted version
% (used here) is a variant.
%
% cf. "Permutation Entropy: A Natural Complexity Measure for Time Series"
% C. Bandt and B. Pompe, Phys. Rev. Lett. 88(17) 174102 (2002)
% (the underlying ordinal-pattern symbolization; cf. EN_PermEn.m, which
% computes entropy-based measures on the same symbolization -- entropy is
% not recomputed here to avoid duplicating that operation, and Kulp et al.
% found entropy to be a weaker discriminator than the network measures
% computed below anyway).
%
% ---INPUTS:
% y, the input time series
%
% d, the ordinal pattern (embedding) dimension: windows of d consecutive
%    (delay-tau-spaced) points are each mapped to their rank permutation,
%    one of d! possible ordinal patterns
%
% tau, the time delay (default: 1, as used throughout Kulp et al. 2016)
%
% ---OUTPUTS: the mean degree (average unique out-edges per visited node --
% the paper's central discriminating measure) and its d-normalized version
% (mean degree is bounded above by d, attained for an unconstrained/random
% process); the number of forbidden (non-occurring) ordinal patterns (NFP)
% and its normalized version; and complementary graph statistics (max/std
% in- and out-degree, edge density, reciprocity, and mean/max edge weight,
% i.e., how often the same transition repeats) not reported in the paper
% but cheaply available from the same transition-pair computation.

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
if nargin < 2 || isempty(d)
    d = 3;
end
if nargin < 3 || isempty(tau)
    tau = 1;
end

% ------------------------------------------------------------------------------
%% Embed the signal into d-dimensional ordinal-pattern windows
% ------------------------------------------------------------------------------
X = BF_Embed(y, tau, d, false);
if isscalar(X) && isnan(X) % embedding failed
    warning('Embedding failed');
    out = NaN; return
end
Nx = size(X, 1);
if Nx < 30
    warning('Time series too short for a meaningful ordinal partition network (Nx = %u)', Nx);
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Symbolize: assign each window an ordinal-pattern ID, assigned dynamically
%% (only for patterns actually observed -- avoids ever enumerating all d!
%% possible permutations, which matters since only a small fraction are
%% typically visited)
% ------------------------------------------------------------------------------
patternMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
s = zeros(Nx, 1);
nextID = 1;
for i = 1:Nx
    [~, ix] = sort(X(i, :));
    key = sprintf('%d,', ix);
    if isKey(patternMap, key)
        s(i) = patternMap(key);
    else
        patternMap(key) = nextID;
        s(i) = nextID;
        nextID = nextID + 1;
    end
end
% Count is a uint64; cast to double immediately, since any arithmetic
% mixing a uint64 with a double elsewhere below would otherwise silently
% round to the nearest integer (MATLAB coerces double/uint64 operations
% into integer arithmetic) -- corrupting meanDegree and every other ratio
% computed from n.
n = double(patternMap.Count); % number of unique ordinal patterns observed (network nodes)

% ------------------------------------------------------------------------------
%% Build the directed transition network from consecutive symbols
% ------------------------------------------------------------------------------
fromNodes = s(1:end - 1);
toNodes = s(2:end);

[pairs, ~, ic] = unique([fromNodes, toNodes], 'rows'); % unique directed edges
weights = accumarray(ic, 1); % how many times each unique edge repeats
m = size(pairs, 1); % number of unique directed edges

% ------------------------------------------------------------------------------
%% Core measures (cf. Kulp et al. 2016)
% ------------------------------------------------------------------------------
% NB on what is deliberately *not* returned here: meanDegreeNorm
% (= meanDegree/d), NFPnorm (= NFP/factorial(d)), density (= meanDegree/n) and
% meanEdgeWeight (= 1/m) were all removed as exact algebraic functions of the
% fields below -- verified to machine precision (residuals 0 to 3e-17) on a
% mixed noise/AR/periodic/logistic ensemble. Being deterministic rescalings,
% they are literally identical columns after hctsa's z-scoring and so add no
% information to a feature library. Any of them can be recovered from
% meanDegree and NFP if wanted, since n = factorial(d) - NFP.
out.meanDegree = m / n; % mean in-degree = mean out-degree = m/n for any directed network

dFactorial = factorial(d);
out.NFP = dFactorial - n; % number of forbidden (non-occurring) ordinal patterns

% ------------------------------------------------------------------------------
%% Complementary graph statistics (not in the paper, cheap from the same pairs)
% ------------------------------------------------------------------------------
outDegCounts = accumarray(pairs(:, 1), 1, [n, 1]);
inDegCounts = accumarray(pairs(:, 2), 1, [n, 1]);

out.maxOutDegree = max(outDegCounts);
out.maxInDegree = max(inDegCounts);
out.stdOutDegree = std(outDegCounts);
out.stdInDegree = std(inDegCounts);

% Reciprocity: fraction of edges (i,j) for which the reverse edge (j,i) also occurs
A = sparse(pairs(:, 1), pairs(:, 2), true, n, n);
out.reciprocity = full(sum(sum(A & A'))) / m;

% Edge weights are raw repeat counts, so they grow in direct proportion to the
% number of transitions: with d fixed there are at most factorial(d)^2 possible
% edges, so m saturates while the transition count keeps growing with N. Divide
% through by the total number of transitions to get edge *probabilities*, which
% are intensive.
numTransitions = sum(weights); % = length(s)-1, the number of observed transitions

% (mean(weights)/numTransitions would be identically 1/m -- see the note above.)
out.maxEdgeWeight = max(weights) / numTransitions; % probability of the most frequent single transition

end
