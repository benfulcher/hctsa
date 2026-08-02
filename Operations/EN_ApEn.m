function out = EN_ApEn(y, mnom, rth)
% EN_ApEn   Approximate Entropy of a time series
%
% ApEn(m,r).
%
% cf. S. M. Pincus, "Approximate entropy as a measure of system complexity",
% P. Natl. Acad. Sci. USA, 88(6) 2297 (1991)
%
% For more information, cf. http://physionet.org/physiotools/ApEn/
%
% ---INPUTS:
% y, the input time series
% mnom, the embedding dimension
% rth, the threshold for judging closeness/similarity
%
% ---NOTES:
% I have no record of where this was code was derived from :-/

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

% Check inputs, set defaults:
if nargin < 2 || isempty(mnom)
	mnom = 1; % m = 1 (default)
end

if nargin < 3 || isempty(rth)
	rth = 0.2; % r = 0.2 (default)
end

% -------------------------------------------------------------------------------

r = rth * std(y); % threshold of similarity
N = length(y); % length of time series
phi = zeros(2, 1); % phi(1)=phi_m, phi(2)=phi_{m+1}

for k = 1:2
	m = mnom + k - 1; % pattern length
	C = zeros(N - m + 1, 1);

	% Form vector sequences x from the time series y: x(i,:) = y(i:i+m-1).
	% Built via one vectorized indexing operation instead of an N-m+1
	% iteration loop:
	idx = (1:N - m + 1)' + (0:m - 1);
	x = y(idx);

	for i = 1:N - m + 1
		% m - m(i,:)-style implicit broadcasting subtracts the row x(i,:)
		% from every row of x, giving the same result as explicitly building
		% ax (formerly done via an inner for-loop over j=1:m per i -- an
		% O(N*m) rebuild on every one of the N-m+1 outer iterations) without
		% that per-iteration loop. The outer loop over i is kept (rather than
		% vectorizing across all i at once) to avoid an O(N^2*m) intermediate
		% array that could be excessive memory for long time series:
		d = abs(x - x(i, :));
		if m > 1 % Takes maximum distance
			d = max(d, [], 2)';
		end
		dr = (d <= r);
		C(i) = sum(dr) / (N - m + 1); % Number of x(j) within r of x(i)
	end
	phi(k) = mean(log(C));
end
out = phi(1) - phi(2);

end
