function out = EN_BubbleEn(y, m, tau)
% EN_BubbleEn   Bubble entropy of a time series.
%
% Manis et al.'s bubble entropy, an ordinal entropy that depends very little on
% its embedding dimension. The series is cut into overlapping runs of m values
% spaced tau samples apart. For each run, the number of swaps a bubble sort needs
% to put it in order (equivalently, the number of pairs in which an earlier value
% exceeds a later one, from 0 to m(m-1)/2) is counted. The Renyi entropy of order
% 2 of the distribution of this swap count, H_m = -ln(sum_k p_k^2), is computed
% for runs of m and of m+1 values. The bubble entropy is the increase in entropy
% on going from m to m+1 values, H_(m+1) - H_m, divided by ln((m+1)/(m-1)) to
% normalize for the dimension. Low values indicate series whose runs have
% predictable orderings. For white noise, every ordering of a run is equally
% likely, and the long-series value follows from the distribution of the number
% of pairs in the wrong order: about 0.64 for m = 5 and 0.69 for m = 10 (rising
% slowly toward 0.75 as m grows).
%
% ---INPUTS:
% y, the input time series
% m, the embedding dimension, at least 2 (default: 10)
% tau, the time delay for the embedding (default: 1); can also be 'ac' (first
%    zero-crossing of the autocorrelation function), 'ac1e' (the floor of its first
%    1/e crossing), or 'mi' (the smaller of the first minimum of the Kraskov
%    automutual information and the 'ac1e' delay), as in BF_Embed (see BF_GetTau)
%
% ---OUTPUTS:
% A structure with a single field:
% bubbleEn, the bubble entropy, (H_(m+1) - H_m)/ln((m+1)/(m-1))
% NaN is returned if the series is too short for runs of m+1 values (fewer than
% 10 runs) or is constant.
%
% ---REFERENCES:
% G. Manis, M. D. Aktaruzzaman and R. Sassi, "Bubble Entropy: An Entropy Almost
% Free of Parameters", IEEE Trans. Biomed. Eng. 64(11), 2711-2718 (2017).
% DOI: 10.1109/TBME.2017.2664105
%
% ---NOTES:
% A swap is counted only when an earlier value is strictly greater than a later
% one, so tied values are never swapped (as in a standard bubble sort).
% The estimate is a small difference between two entropies, so it is noisy when
% the series is short relative to the number of possible swap counts: for series
% of about 1000 samples, embedding dimensions much above 10 give poorly
% reproducible values.
%
% With a delay set by the autocorrelation function ('ac'), the runs of m values
% span (m-1)*tau samples, so for a slowly decorrelating series (a large delay)
% there are few runs and the value is unreliable; it is NaN when the
% autocorrelation function has no zero crossing (as for many random walks).

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
if nargin < 2 || isempty(m), m = 10; end
if nargin < 3 || isempty(tau), tau = 1; end

out.bubbleEn = NaN;

y = y(:);
if ~(std(y) > 0)
	return % constant (or non-finite) series
end
if ischar(tau) % resolve 'ac' / 'mi' once, so both dimensions share one delay
	p = BF_Embed(y, tau, m, true);
	if isscalar(p) && isnan(p), return; end
	tau = p(1);
end

% ------------------------------------------------------------------------------
%% Renyi-2 entropy of the swap-count distribution at m and m+1:
% ------------------------------------------------------------------------------
H = zeros(1, 2);
for j = 1:2
	mm = m + j - 1;
	x = BF_Embed(y, tau, mm, false);
	if isscalar(x) || size(x, 1) < 10
		return
	end
	numSwaps = zeros(size(x, 1), 1);
	for a = 1:mm - 1
		for b = a + 1:mm
			numSwaps = numSwaps + (x(:, a) > x(:, b)); % one swap per inversion
		end
	end
	p = accumarray(numSwaps + 1, 1) / size(x, 1);
	H(j) = -log(sum(p.^2));
end

out.bubbleEn = (H(2) - H(1)) / log((m + 1) / (m - 1));

end
