function out = EN_RandomizeOrdinal(y, tau)
% EN_RandomizeOrdinal   How quickly ordinal-pattern structure is lost as a series is randomized.
%
% The series is randomized progressively: at each step, one randomly chosen
% sample is overwritten by a value copied from a randomly chosen sample of the
% original series (the 'statdist' scheme of EN_Randomize). Rather than simulating
% this, the function computes exactly (with no random draws) the expected
% distribution of the ordinal patterns of order 3 (the patterns of Bandt and
% Pompe's permutation entropy, PermEn(3,tau)) after t steps, and reports how
% quickly the ordinal structure of the series is lost.
%
% The calculation is exact because, after t steps, each sample has either never
% been overwritten, or it holds a value drawn independently and uniformly from the
% original series. For a window of three samples, the probability that exactly a
% given subset of k of them has never been overwritten is
%       w_k(t) = sum_{j=0}^{3-k} nchoosek(3-k, j) (-1)^j (1 - (k+j)/N)^t.
% Given the subset, the pattern probabilities follow from the kept values and the
% empirical distribution function of the series (ties are treated as in
% EN_PermEn: equal values are ordered by their position). The expected pattern
% distribution after t steps is therefore a fixed combination of eight
% series-specific pattern distributions (one per subset), weighted by w_k(t),
% and costs O(N log N) to compute.
%
% The main output, halftime, is the randomization time t/N (the number of
% replacements per sample) at which the permutation entropy of the expected
% pattern distribution is halfway between its value for the original series
% (t = 0) and its value for the fully randomized series (t -> infinity). A long
% halftime means that much of the ordinal structure survives the replacement of
% single samples (it is carried by the order of pairs of samples); a short
% halftime means that the structure depends on the joint configuration of all
% three samples of a window.
%
% ---INPUTS:
% y, the input time series
% tau, the time delay between the three samples of each window (default: 1); can
%    also be a rule understood by BF_GetTau ('ac', 'ac1e', 'mi'), evaluated on y
%
% ---OUTPUTS:
% A structure with fields:
% halftime, the randomization time t/N at which the permutation entropy of the
%    expected pattern distribution is halfway to its fully randomized value (NaN
%    if the original series has the same pattern entropy as its fully randomized
%    version, so that nothing decays)
% normPermEn0, the normalized PermEn(3,tau) of the original series
% normPermEnRep1, the normalized permutation entropy of the expected pattern
%    distribution when one sample of each window (chosen at random) is replaced by
%    an independent draw from the series
% normPermEnInf, the normalized permutation entropy of the fully randomized series
%    (close to 1 for a series without repeated values)
% NaN (instead of a structure) is returned if the series is too short to embed
% (fewer than 5 windows).
%
% ---NOTES:
% The entropy is that of the expected pattern distribution (not the expected
% entropy of a single randomized series, which is lower by a small finite-sample
% bias), so the result is deterministic and needs no random seed. The 'dyndist'
% and 'permute' schemes of EN_Randomize give the same expected pattern
% distribution to leading order in 1/N ('permute' at twice the rate), so they
% are not offered separately.
%
% ---REFERENCES:
% C. Bandt and B. Pompe, "Permutation Entropy: A Natural Complexity Measure for
% Time Series", Phys. Rev. Lett. 88(17) 174102 (2002).
% DOI: 10.1103/PhysRevLett.88.174102

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

if nargin < 2 || isempty(tau)
	tau = 1;
end

y = y(:);
N = length(y);
if ischar(tau)
	tau = BF_GetTau(y, tau); % the delay is set by the series (not by its ranks)
	if isnan(tau)
		out = NaN; return
	end
end
if N - 2 * tau < 5 % need at least 5 windows [y(i), y(i+tau), y(i+2*tau)]
	warning('Time series too short to embed');
	out = NaN; return
end

% ------------------------------------------------------------------------------
% Pattern distributions for each subset of kept samples
% ------------------------------------------------------------------------------
[S, subsetSize] = SUB_PatternDists(y, tau);

H = @(P) SUB_NormEntropy(P);
out.normPermEn0 = H(S(subsetSize == 3, :));
out.normPermEnRep1 = H(mean(S(subsetSize == 2, :), 1));
out.normPermEnInf = H(S(subsetSize == 0, :));

% ------------------------------------------------------------------------------
% Half time: first crossing of the midpoint entropy, on a grid then by bisection
% ------------------------------------------------------------------------------
H0 = out.normPermEn0;
Hinf = out.normPermEnInf;
if abs(Hinf - H0) < 1e-12
	out.halftime = NaN;
	return
end
target = 0.5 * (H0 + Hinf);
sgn = sign(Hinf - H0);
g = @(tN) sgn * (H(SUB_Weights(N, tN * N, subsetSize) * S) - target);

gridStep = 0.01; tMax = 10; % in units of t/N
tGrid = (0:gridStep:tMax)';
gGrid = g(tGrid);
j = find(gGrid >= 0, 1, 'first');
if isempty(j)
	out.halftime = NaN;
	return
end
lo = tGrid(j - 1); hi = tGrid(j);
for k = 1:60
	mid = 0.5 * (lo + hi);
	if g(mid) >= 0
		hi = mid;
	else
		lo = mid;
	end
end
out.halftime = 0.5 * (lo + hi);

end
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
function [S, subsetSize] = SUB_PatternDists(y, tau)
	% S(k, :): expected ordinal-pattern distribution (averaged over windows) when
	% the samples in subset k of each window are kept and the others are
	% independent draws from the empirical distribution of y
	N = length(y);
	[~, ~, r] = unique(y); % ranks 1..U of the distinct values
	U = max(r);
	f = accumarray(r, 1, [U, 1]) / N; % probability of each distinct value
	F = cumsum(f); % P(X <= k)
	F2 = cumsum(f.^2);
	F0 = [0; F]; F20 = [0; F2];
	cdf = @(k) F0(k + 1); % P(X <= k) for k = 0..U
	cdf2 = @(k) F20(k + 1);
	pm = @(k) f(k);
	W = BF_Embed(r, tau, 3, false); % windows of ranks (same windows as BF_Embed(y, tau, 3))
	Nx = size(W, 1);

	subsets = {[], 1, 2, 3, [1, 2], [1, 3], [2, 3], [1, 2, 3]};
	subsetSize = cellfun(@numel, subsets);
	S = zeros(8, 6);
	for k = 1:8
		K = subsets{k};
		free = setdiff(1:3, K);
		switch numel(K)
		case 3
			S(k, :) = accumarray(SUB_Pattern(W(:, 1), W(:, 2), W(:, 3)), 1, [6, 1])';
		case 2
			a = K(1); b = K(2); c = free;
			lo = min(W(:, a), W(:, b)); hi = max(W(:, a), W(:, b)); same = (lo == hi);
			xv = {lo - 0.5, lo, lo + 0.5, hi, hi + 0.5};
			pr = {cdf(lo - 1), pm(lo), (~same) .* (cdf(hi - 1) - cdf(lo)), (~same) .* pm(hi), 1 - cdf(hi)};
			for q = 1:5
				Z = cell(1, 3); Z{a} = W(:, a); Z{b} = W(:, b); Z{c} = xv{q} + zeros(Nx, 1);
				S(k, :) = S(k, :) + accumarray(SUB_Pattern(Z{:}), pr{q}, [6, 1])';
			end
		case 1
			a = K; b = free(1); c = free(2); v = W(:, a);
			pL = cdf(v - 1); pE = pm(v); pG = 1 - cdf(v);
			sL = cdf2(v - 1); sG = F2(end) - cdf2(v);
			regionP = {pL, pE, pG}; regionRep = {v - 0.5, v, v + 0.5};
			for R1 = 1:3
				for R2 = 1:3
					if R1 == R2 && R1 ~= 2
						% both draws in the same strict region: order between them
						if R1 == 1
							p = pL; s = sL; x1 = v - 0.75; x2 = v - 0.5;
						else
							p = pG; s = sG; x1 = v + 0.5; x2 = v + 0.75;
						end
						half = (p.^2 - s) / 2;
						cases = {x1, x2, half; x2, x1, half; x1, x1, s};
						for q = 1:3
							Z = cell(1, 3); Z{a} = v; Z{b} = cases{q, 1}; Z{c} = cases{q, 2};
							S(k, :) = S(k, :) + accumarray(SUB_Pattern(Z{:}), cases{q, 3}, [6, 1])';
						end
					else
						Z = cell(1, 3); Z{a} = v; Z{b} = regionRep{R1}; Z{c} = regionRep{R2};
						S(k, :) = S(k, :) + accumarray(SUB_Pattern(Z{:}), regionP{R1} .* regionP{R2}, [6, 1])';
					end
				end
			end
		case 0
			% three independent draws: the same for every window
			s2 = sum(f.^2); s3 = sum(f.^3);
			pDistinct = 1 - 3 * (s2 - s3) - s3;
			q = pDistinct / 6 * ones(1, 6);
			below = sum(f.^2 .* [0; F(1:end-1)]); % two equal draws, the third below them
			above = sum(f.^2 .* (1 - F)); % ... or above them
			pairs = [1, 2; 1, 3; 2, 3];
			for p = 1:3
				o = setdiff(1:3, pairs(p, :));
				z = zeros(1, 3); z(o) = -1; q(SUB_Pattern(z(1), z(2), z(3))) = q(SUB_Pattern(z(1), z(2), z(3))) + below;
				z = zeros(1, 3); z(o) = 1; q(SUB_Pattern(z(1), z(2), z(3))) = q(SUB_Pattern(z(1), z(2), z(3))) + above;
			end
			q(SUB_Pattern(0, 0, 0)) = q(SUB_Pattern(0, 0, 0)) + s3;
			S(k, :) = q * Nx;
		end
	end
	S = S / Nx;
end
% ------------------------------------------------------------------------------
function idx = SUB_Pattern(z1, z2, z3)
	% Ordinal pattern index (as BF_OrdinalPatternRank: stable sort, so equal values
	% are ordered by position) from the three pairwise precedence relations
	persistent LUT
	if isempty(LUT)
		LUT = zeros(8, 1);
		P = perms(1:3);
		for i = 1:size(P, 1)
			p = P(i, :);
			code = 4 * (p(1) <= p(2)) + 2 * (p(1) <= p(3)) + (p(2) <= p(3));
			LUT(code + 1) = BF_OrdinalPatternRank(p);
		end
	end
	code = 4 * (z1 <= z2) + 2 * (z1 <= z3) + (z2 <= z3);
	idx = LUT(code + 1);
end
% ------------------------------------------------------------------------------
function w = SUB_Weights(N, t, subsetSize)
	% w(i, k): probability that, after t(i) steps, exactly the samples of subset k
	% of a window have never been overwritten
	t = t(:);
	w = zeros(length(t), length(subsetSize));
	for k = 1:length(subsetSize)
		s = subsetSize(k);
		for j = 0:3 - s
			w(:, k) = w(:, k) + nchoosek(3 - s, j) * (-1)^j * (1 - (s + j) / N).^t;
		end
	end
end
% ------------------------------------------------------------------------------
function h = SUB_NormEntropy(P)
	% normalized Shannon entropy (base 6) of each row of P
	Pl = P; Pl(Pl <= 0) = 1;
	h = -sum(P .* log2(Pl), 2) / log2(6);
end
