function x = BF_Random(n, seed, kind)
% BF_Random  Portable pseudo-random numbers: the same stream in MATLAB and Python.
%
% A small implementation of L'Ecuyer's combined multiple recursive generator MRG32k3a
% (L'Ecuyer 1999). It uses double-precision arithmetic only (every product stays below
% 2^53, so there is no integer overflow), which makes it a few lines in any language
% and gives bit-identical uniform streams everywhere. It leaves MATLAB's global random
% stream untouched, so using it neither disturbs nor depends on rng.
%
% Used where a feature needs random numbers that are reproducible across MATLAB
% releases and across implementations (MATLAB's randn and randperm are not).
%
% ---INPUTS:
% n, the number of values to return
% seed, an integer in [0, 4e9) that sets the stream (default: 0). All six state words
%       are set to 12345 + seed (12345 in all six is the reference generator's default
%       state) and the first 8 draws are discarded, so that streams from neighboring
%       seeds are unrelated. A vector of six state words [s10 s11 s12 s20 s21 s22]
%       (the first three in [0, 4294967087), the last three in [0, 4294944443),
%       neither triple all zero) instead starts the raw generator with no draws
%       discarded; the state [12345 12345 12345 12345 12345 12345] reproduces
%       L'Ecuyer's reference implementation exactly.
% kind, what to return:
%       'uniform': n uniform numbers in the open interval (0,1) (default)
%       'normal': n standard normal numbers by the Box-Muller transform: the
%                 uniforms (u1, u2) = (2k-1, 2k) give the normal numbers
%                 sqrt(-2*log(u1))*cos(2*pi*u2) and sqrt(-2*log(u1))*sin(2*pi*u2)
%                 as values 2k-1 and 2k
%       'perm': a random permutation of 1..n: the ranks of n uniform numbers, i.e.
%               the indices that sort them ([~, x] = sort(u); ties, which have a
%               probability of about n^2/2^33, go to the lower index)
%
% ---OUTPUTS:
% x, a column vector (n-by-1)
%
% ---NOTES:
% With state s1 (s10, s11, s12) and s2 (s20, s21, s22), each draw computes
%       p1 = (1403580*s11 - 810728*s10) mod 4294967087, s1 <- (s11, s12, p1)
%       p2 = (527612*s22 - 1370589*s20) mod 4294944443, s2 <- (s21, s22, p2)
% and returns (p1 - p2 [+ 4294967087 if p1 <= p2]) * 2.328306549295727688e-10.
% Uniform draws agree exactly between languages; normal draws agree to the
% rounding of the platform's log, cos and sin (about 1e-16).
%
% Speed: the recurrence is linear, so a long stream is cut into blocks of 2^6 draws
% whose starting states are obtained by jump-ahead (powers of the 3-by-3 transition
% matrices modulo m1 and m2, in exact uint64 arithmetic) and the blocks are then
% advanced together as columns. This returns exactly the same numbers as the plain
% one-draw-at-a-time loop (used for short streams), only faster for long ones.
%
% ---REFERENCES:
% P. L'Ecuyer, "Good parameters and implementations for combined multiple recursive
% random number generators", Operations Research 47(1), 159 (1999).
% DOI: 10.1287/opre.47.1.159
% P. L'Ecuyer, R. Simard, E. J. Chen and W. D. Kelton, "An object-oriented
% random-number package with many long streams and substreams", Operations Research
% 50(6), 1073 (2002). DOI: 10.1287/opre.50.6.1073.358

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

if nargin < 2 || isempty(seed)
	seed = 0;
end
if nargin < 3 || isempty(kind)
	kind = 'uniform';
end

if isscalar(seed)
	s0 = (12345 + seed) * ones(1, 6); % the six state words
	numSkip = 8; % decorrelate neighboring seeds
else
	s0 = seed(:)'; % raw state
	numSkip = 0;
end

switch kind
	case 'uniform'
		x = nextUniforms(s0, numSkip, n);
	case 'normal'
		m = ceil(n / 2);
		u = nextUniforms(s0, numSkip, 2 * m);
		r = sqrt(-2 * log(u(1:2:end)));
		theta = 2 * pi * u(2:2:end);
		x = reshape([r .* cos(theta), r .* sin(theta)]', [], 1);
		x = x(1:n);
	case 'perm'
		[~, x] = sort(nextUniforms(s0, numSkip, n));
	otherwise
		error('Unknown kind ''%s''', kind);
end

end

% ------------------------------------------------------------------------------
function u = nextUniforms(s0, numSkip, k)
% k uniforms after discarding numSkip draws, starting from the state s0 (six words)
m1 = 4294967087;
m2 = 4294944443;
total = k + numSkip;
if total <= 3000
	p1 = zeros(total, 1);
	p2 = zeros(total, 1);
	a1 = s0(1); b1 = s0(2); c1 = s0(3);
	a2 = s0(4); b2 = s0(5); c2 = s0(6);
	for ii = 1:total
		d1 = mod(1403580 * b1 - 810728 * a1, m1);
		a1 = b1; b1 = c1; c1 = d1;
		d2 = mod(527612 * c2 - 1370589 * a2, m2);
		a2 = b2; b2 = c2; c2 = d2;
		p1(ii) = d1;
		p2(ii) = d2;
	end
else
	% blocks of L = 2^p draws, B of them side by side
	p = 6;
	L = 2^p;
	B = ceil(total / L);
	J1 = jumpMatrices(1);
	J2 = jumpMatrices(2);
	X1 = uint64(s0(1:3)');
	X2 = uint64(s0(4:6)');
	q = 0;
	while size(X1, 2) < B % block j starts at A^(jL) * s0, built up by doubling
		X1 = [X1, applyMod(J1{p + q + 1}, X1, uint64(m1))]; %#ok<AGROW>
		X2 = [X2, applyMod(J2{p + q + 1}, X2, uint64(m2))]; %#ok<AGROW>
		q = q + 1;
	end
	a1 = double(X1(1, 1:B))'; b1 = double(X1(2, 1:B))'; c1 = double(X1(3, 1:B))';
	a2 = double(X2(1, 1:B))'; b2 = double(X2(2, 1:B))'; c2 = double(X2(3, 1:B))';
	P1 = zeros(B, L);
	P2 = zeros(B, L);
	for ii = 1:L
		% (x - m*floor(x/m), corrected by one m, is exact here and much faster than mod)
		d1 = 1403580 * b1 - 810728 * a1;
		d1 = d1 - m1 * floor(d1 / m1);
		d1 = d1 + m1 * (d1 < 0) - m1 * (d1 >= m1);
		a1 = b1; b1 = c1; c1 = d1;
		d2 = 527612 * c2 - 1370589 * a2;
		d2 = d2 - m2 * floor(d2 / m2);
		d2 = d2 + m2 * (d2 < 0) - m2 * (d2 >= m2);
		a2 = b2; b2 = c2; c2 = d2;
		P1(:, ii) = d1;
		P2(:, ii) = d2;
	end
	p1 = P1.'; % column j holds draws (j-1)*L+1 ... j*L
	p2 = P2.';
	p1 = p1(1:total)';
	p2 = p2(1:total)';
end
p1 = p1(numSkip + 1:end);
p2 = p2(numSkip + 1:end);
u = (p1 - p2 + m1 * (p1 <= p2)) * 2.328306549295727688e-10;
u = reshape(u, [], 1);
end

% ------------------------------------------------------------------------------
function J = jumpMatrices(comp)
% J{q+1} = (transition matrix of component comp)^(2^q) modulo its modulus, exact
persistent cache
if isempty(cache)
	cache = cell(2, 1);
end
if isempty(cache{comp})
	if comp == 1
		m = uint64(4294967087);
		A = uint64([0 1 0; 0 0 1; 4294967087 - 810728, 1403580, 0]);
	else
		m = uint64(4294944443);
		A = uint64([0 1 0; 0 0 1; 4294944443 - 1370589, 0, 527612]);
	end
	J = cell(1, 40);
	J{1} = A;
	for q = 2:40
		J{q} = applyMod(J{q - 1}, J{q - 1}, m);
	end
	cache{comp} = J;
end
J = cache{comp};
end

function Y = applyMod(M, X, m)
% M*X modulo m for uint64 M (3-by-3) and X (3-by-k), entries below m (< 2^32), so every
% product is exact in uint64
Y = zeros(3, size(X, 2), 'uint64');
for r = 1:3
	Y(r, :) = rem(rem(M(r, 1) * X(1, :), m) + rem(M(r, 2) * X(2, :), m) + rem(M(r, 3) * X(3, :), m), m);
end
end

