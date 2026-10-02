function out = EN_OrdinalPatterns(y, m, tau)
% EN_OrdinalPatterns   Amplitude-weighted entropy, forbidden patterns and time asymmetry of ordinal patterns.
%
% Cuts the series into overlapping runs of m values spaced tau samples apart and
% replaces each run by its ordinal pattern (the order its m values would take if
% sorted, one of m! possible patterns), as in EN_PermEn. Three properties of the
% resulting pattern distribution are returned that EN_PermEn does not give:
% (i) the weighted permutation entropy, in which each run counts in proportion to
% its variance, so patterns formed by large fluctuations matter more than those
% formed by small ones; (ii) the fraction of the m! patterns that never occur
% ('forbidden' patterns, which deterministic dynamics have and noise, given
% enough data, does not); and (iii) a time asymmetry, the distance between the
% pattern distribution of the series and that of the same series read backward.
%
% ---INPUTS:
% y, the input time series
% m, the embedding dimension (the order of the patterns; default: 4)
% tau, the time delay for the embedding (default: 1); can also be 'ac' (first
%    zero-crossing of the autocorrelation function) or 'mi' (first minimum of
%    the automutual information), as in BF_Embed
%
% ---OUTPUTS:
% A structure with fields:
% normWPE, the weighted permutation entropy, normalized by log2(m!), from 0 to 1;
%    each pattern's probability is the sum of the variances of the runs that have
%    that pattern, divided by the sum of the variances of all runs
% forbidFrac, the fraction of the m! possible ordinal patterns that occur in no
%    run at all, from 0 to 1
% ordAsym, the time asymmetry: the total variation distance,
%    (1/2) sum_j |p_j - q_j|, between the pattern distribution p of the series
%    and the pattern distribution q of the time-reversed series, from 0 (equal
%    distributions) to 1
% All fields are NaN if the series is constant or too short to embed (fewer than
% 5 embedding vectors).
%
% ---REFERENCES:
% B. Fadlallah, B. Chen, A. Keil and J. Principe, "Weighted-permutation
% entropy: A complexity measure for time series incorporating amplitude
% information", Phys. Rev. E 87, 022911 (2013). DOI: 10.1103/PhysRevE.87.022911
%
% J. M. Amigo, S. Zambrano and M. A. F. Sanjuan, "True and false forbidden
% patterns in deterministic and random dynamics", Europhys. Lett. 79, 50001
% (2007). DOI: 10.1209/0295-5075/79/50001
%
% M. Zanin, A. Rodriguez-Gonzalez, E. Menasalvas Ruiz and D. Papo, "Assessing
% time series reversibility through permutation patterns", Entropy 20(9), 665
% (2018). DOI: 10.3390/e20090665
%
% ---NOTES:
% The time-reversed series has the same runs as the original read backward, so
% its pattern distribution is computed from the same runs without re-embedding.
% ordAsym is a distance between two empirical distributions and so is above zero
% even for a time-reversible series; the amount is about sqrt(m!/Nx) for Nx runs.
% forbidFrac depends on the length of the series: it is only informative when the
% number of runs is much larger than m!.

% ------------------------------------------------------------------------------
% Copyright (C) 2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
%% Check inputs and set defaults:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(m), m = 4; end
if nargin < 3 || isempty(tau), tau = 1; end

out.normWPE = NaN;
out.forbidFrac = NaN;
out.ordAsym = NaN;

y = y(:);
if ~(std(y) > 0)
	return % constant (or non-finite) series
end

% ------------------------------------------------------------------------------
%% Embed the signal and find the ordinal pattern of each run (as in EN_PermEn):
% ------------------------------------------------------------------------------
x = BF_Embed(y, tau, m, false);
Nx = size(x, 1); % number of embedding vectors produced
if isscalar(x) || Nx < 5 % too short for the requested tau/m
	return
end

numPerms = factorial(m);
permIdx = BF_OrdinalPatternRank(x);
countPerms = accumarray(permIdx, 1, [numPerms, 1]);

% ------------------------------------------------------------------------------
%% Weighted permutation entropy (each run weighted by its variance):
% ------------------------------------------------------------------------------
w = var(x, 1, 2); % variance of the m values in each run
if sum(w) > 0
	pw = accumarray(permIdx, w, [numPerms, 1]) / sum(w);
	pw = pw(pw > 0);
	out.normWPE = -sum(pw .* log2(pw)) / log2(numPerms);
end

% ------------------------------------------------------------------------------
%% Fraction of forbidden patterns:
% ------------------------------------------------------------------------------
out.forbidFrac = mean(countPerms == 0);

% ------------------------------------------------------------------------------
%% Ordinal time asymmetry (the reversed series has each run read backward):
% ------------------------------------------------------------------------------
permIdxRev = BF_OrdinalPatternRank(x(:, end:-1:1));
countPermsRev = accumarray(permIdxRev, 1, [numPerms, 1]);
out.ordAsym = 0.5 * sum(abs(countPerms - countPermsRev)) / Nx;

end
