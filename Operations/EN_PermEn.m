function out = EN_PermEn(y, m, tau)
% EN_PermEn   Permutation entropy of a time series.
%
% Bandt and Pompe's permutation entropy. The series is cut into overlapping runs
% of m values, spaced tau samples apart, and each run is replaced by its ordinal
% pattern (the order its m values would take if sorted, one of m! possible
% patterns). The permutation entropy is the Shannon entropy of how often each
% pattern occurs; it depends only on the order of the values, not their size.
% Also returns the version normalized by log2(m!), and an adapted implementation
% by Bruce Land and Damian Elias.
%
% ---INPUTS:
% y, the input time series
% m, the embedding dimension (the order of the permutation entropy; default: 2)
% tau, the time delay for the embedding (default: 1); can also be 'ac' (first
%    zero-crossing of the autocorrelation function) or 'mi' (first minimum of
%    the automutual information), as in BF_Embed
%
% ---OUTPUTS:
% A structure with fields:
% permEn, the permutation entropy in bits, -sum(p.*log2(p)) over the ordinal-
%    pattern probabilities p
% normPermEn, permEn normalized by log2(m!), from 0 to 1
% permEnLE, the Land-Elias version: the entropy in nats, with probabilities below
%    1/(number of embedding vectors) raised to that floor, divided by (m - 1)
% NaN (instead of a structure) is returned if the series is too short to embed
% (fewer than 5 embedding vectors).
%
% ---REFERENCES:
% C. Bandt and B. Pompe, "Permutation Entropy: A Natural Complexity Measure for
% Time Series", Phys. Rev. Lett. 88(17) 174102 (2002).
%
% ---NOTES:
% The Land-Elias version is adapted from
% http://people.ece.cornell.edu/land/PROJECTS/Complexity/ (logisticPE.m).

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

% -------------------------------------------------------------------------------
% Check inputs and set defaults:
% -------------------------------------------------------------------------------

if nargin < 2 || isempty(m)
	m = 2; % order 2
end

if nargin < 3 || isempty(tau)
	tau = 1;
end

% ------------------------------------------------------------------------------
% Embed the signal:
% ------------------------------------------------------------------------------
x = BF_Embed(y, tau, m, false);
Nx = size(x, 1); % number of embedding vectors produced
if Nx < 5 % need at least 5 embedding vectors to actually do a computation
	% Data-dependent (series too short for the requested tau/m), not a code error.
	warning('Time series too short to embed');
	out = NaN; return
end
numPerms = factorial(m);
permIdx = BF_OrdinalPatternRank(x); % index in 1:m! for each embedding vector
countPerms = accumarray(permIdx, 1, [numPerms, 1]);

% ------------------------------------------------------------------------------
% Convert counts to probabilities
p = countPerms / Nx; % ((Nx-(m-1))*tau);
p_0 = p(p > 0); % makes log(0) = 0
out.permEn = -sum(p_0 .* log2(p_0));

% Normalized permutation entropy (more comparable across m?)
mFact = factorial(m);
out.normPermEn = out.permEn / log2(mFact);

% ------------------------------------------------------------------------------
% Adapted implementation by Bruce Land and Damian Elias:
% cf.:
% http://people.ece.cornell.edu/land/PROJECTS/Complexity/
% http://people.ece.cornell.edu/land/PROJECTS/Complexity/logisticPE.m

% Not clear to me why you would make log(0) = log(1/Nx); (the minimum)
% rather than exclude it from the sum, as is done here:
p_LE = max(1 / Nx, p);
out.permEnLE = -sum(p_LE .* log(p_LE)) / (m - 1);

end
