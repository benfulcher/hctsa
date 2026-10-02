function out = EN_FuzzyEn(y, M, r, n)
% EN_FuzzyEn   Fuzzy entropy of a time series.
%
% Chen et al.'s fuzzy entropy, a smooth relative of sample entropy. The series is
% cut into overlapping runs of m consecutive values, and the mean of each run is
% subtracted from it, so runs are compared by shape and not by level. Two runs
% are not simply 'matching' or 'not matching', as in sample entropy (EN_SampEn):
% they are given a similarity exp(-(d/r)^n), where d is the largest absolute
% difference between corresponding values of the two (baseline-removed) runs.
% phi_m is the mean similarity over all pairs of distinct runs of length m, and
% the fuzzy entropy at dimension m is FuzzyEn(m) = ln(phi_m) - ln(phi_(m+1)).
% Low values indicate regular, predictable series; high values irregular ones.
% The smooth similarity makes the measure continuous in r and defined for short
% series for which sample entropy would find no matches.
%
% All runs of length 1, ..., M+1 are taken from the same N-M starting points, so
% that successive dimensions are compared on the same footing. The distances are
% computed in blocks, so memory use does not grow with the square of the series
% length (the run time does: it is O(N^2 M)).
%
% ---INPUTS:
% y, the input time series
% M, the largest embedding dimension: FuzzyEn is returned for m = 1, ..., M
%    (default: 2)
% r, the width of the similarity function, as a fraction of the standard
%    deviation of y (default: 0.2). The width in the units of y is r*std(y), so
%    the measure is unchanged by any rescaling of y.
% n, the exponent of the similarity function exp(-(d/r)^n) (default: 2; larger
%    values make the similarity closer to a hard threshold)
%
% ---OUTPUTS:
% A structure with fields fuzzyEn1, fuzzyEn2, ..., fuzzyEnM, the fuzzy entropy
% at each embedding dimension m (in nats). At m = 1 the run mean removed is the
% value itself, so phi_1 = 1 and fuzzyEn1 = -ln(phi_2). All fields are NaN if the
% series is constant, has fewer than M+3 points, or has no pair of runs with a
% nonzero similarity.
%
% ---REFERENCES:
% W. Chen, Z. Wang, H. Xie and W. Yu, "Characterization of surface EMG signal
% based on fuzzy entropy", IEEE Trans. Neural Syst. Rehabil. Eng. 15(2), 266-272
% (2007). DOI: 10.1109/TNSRE.2007.897025
%
% W. Chen, J. Zhuang, W. Yu and Z. Wang, "Measuring complexity using FuzzyEn,
% ApEn, and SampEn", Med. Eng. Phys. 31(1), 61-68 (2009).
% DOI: 10.1016/j.medengphy.2008.04.005
%
% ---NOTES:
% To get the fuzzy entropy of the increments of a series, give diff(y) as the
% input; r is relative to the standard deviation of whatever series is given.

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
if nargin < 2 || isempty(M), M = 2; end
if nargin < 3 || isempty(r), r = 0.2; end
if nargin < 4 || isempty(n), n = 2; end

y = y(:);
N = length(y);

% Fixed set of output fields (NaN until computed):
fieldNames = arrayfun(@(m) sprintf('fuzzyEn%u', m), 1:M, 'UniformOutput', false);
for i = 1:M
	out.(fieldNames{i}) = NaN;
end

sd = std(y);
Nv = N - M; % number of starting points shared by every embedding dimension
if ~isfinite(sd) || sd == 0 || Nv < 3
	return
end
width = r * sd;

% ------------------------------------------------------------------------------
%% Mean similarity of pairs of distinct baseline-removed runs, at m = 1, ..., M+1:
% ------------------------------------------------------------------------------
phi = zeros(M + 1, 1);
blockSize = max(1, floor(2e6 / Nv)); % cap the size of the distance block
for m = 1:M + 1
	Z = zeros(Nv, m);
	for k = 1:m
		Z(:, k) = y(k:k + Nv - 1);
	end
	Z = Z - mean(Z, 2); % remove each run's own mean (local baseline)

	total = 0;
	for i0 = 1:blockSize:Nv
		ii = i0:min(i0 + blockSize - 1, Nv);
		D = zeros(length(ii), Nv);
		for k = 1:m
			D = max(D, abs(Z(ii, k) - Z(:, k)')); % Chebyshev distance
		end
		total = total + sum(exp(-(D / width).^n), 'all') - length(ii); % drop self-similarity (=1)
	end
	phi(m) = total / (Nv * (Nv - 1));
end

% ------------------------------------------------------------------------------
%% Fuzzy entropy at each dimension:
% ------------------------------------------------------------------------------
if any(phi <= 0)
	return
end
for m = 1:M
	out.(fieldNames{m}) = log(phi(m)) - log(phi(m + 1));
end

end
