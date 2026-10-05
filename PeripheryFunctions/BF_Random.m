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
%       'perm': a random permutation of 1..n by a Fisher-Yates shuffle: for
%               i = n down to 2, element i is swapped with element floor(u*i)+1,
%               where u is the next uniform
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

m1 = 4294967087;
m2 = 4294944443;
if isscalar(seed)
	s = (12345 + seed) * ones(1, 6); % the six state words
	numSkip = 8; % decorrelate neighboring seeds
else
	s = seed(:)'; % raw state
	numSkip = 0;
end

switch kind
	case 'uniform'
		x = nextUniforms(n);
	case 'normal'
		m = ceil(n / 2);
		u = nextUniforms(2 * m);
		r = sqrt(-2 * log(u(1:2:end)));
		theta = 2 * pi * u(2:2:end);
		x = reshape([r .* cos(theta), r .* sin(theta)]', [], 1);
		x = x(1:n);
	case 'perm'
		x = (1:n)';
		u = nextUniforms(max(n - 1, 0));
		for i = n:-1:2
			j = floor(u(n - i + 1) * i) + 1; % uniform on 1..i
			t = x(i); x(i) = x(j); x(j) = t;
		end
	otherwise
		error('Unknown kind ''%s''', kind);
end

	function u = nextUniforms(k)
		u = zeros(k, 1);
		for ii = 1:(k + numSkip)
			p1 = mod(1403580 * s(2) - 810728 * s(1), m1);
			s(1:3) = [s(2), s(3), p1];
			p2 = mod(527612 * s(6) - 1370589 * s(4), m2);
			s(4:6) = [s(5), s(6), p2];
			if ii > numSkip
				if p1 > p2
					u(ii - numSkip) = (p1 - p2) * 2.328306549295727688e-10;
				else
					u(ii - numSkip) = (p1 - p2 + m1) * 2.328306549295727688e-10;
				end
			end
		end
	end

end
