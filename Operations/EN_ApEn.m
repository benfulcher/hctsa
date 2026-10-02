function out = EN_ApEn(y, mnom, rth, tau)
% EN_ApEn   Approximate entropy of a time series.
%
% Pincus's approximate entropy, ApEn(m,r). Every run of m consecutive values is
% compared with every other run; two runs are similar when no pair of
% corresponding values differs by more than r = rth*std(y) (maximum-norm
% distance). Phi_m is the mean over runs of the log-proportion of runs similar to
% each run (a run counts as similar to itself), and the output is
% Phi_m - Phi_{m+1}. Low values indicate regular, predictable series; high
% values irregular ones. Because self-matches are counted, it is biased towards
% low values for short series.
%
% ---INPUTS:
% y, the input time series
% mnom, the embedding dimension m (default: 1)
% rth, the similarity threshold as a fraction of the standard deviation of y,
%      r = rth*std(y) (default: 0.2)
% tau, the time delay between pattern elements (default: 1), or 'ac1e' or 'mi' for an
%      adaptive delay (see BF_GetTau). A delay set by the series' own timescale stops ApEn
%      from mostly measuring smoothness when a process is oversampled, without shortening
%      the series as decimation would.
%
% ---OUTPUTS:
% a scalar: ApEn(m,r) = Phi_m - Phi_{m+1}.
%
% ---REFERENCES:
% S. M. Pincus, "Approximate entropy as a measure of system complexity",
% P. Natl. Acad. Sci. USA, 88(6) 2297 (1991).
%
% ---NOTES:
% For more information, cf. http://physionet.org/physiotools/ApEn/
% I have no record of where this code was derived from :-/

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

% Check inputs, set defaults:
if nargin < 2 || isempty(mnom)
	mnom = 1; % m = 1 (default)
end

if nargin < 3 || isempty(rth)
	rth = 0.2; % r = 0.2 (default)
end

if nargin < 4 || isempty(tau)
	tau = 1; % consecutive samples (default)
elseif ischar(tau)
	tau = BF_GetTau(y, tau);
	if isnan(tau)
		out = NaN; return
	end
end

% -------------------------------------------------------------------------------

r = rth * std(y); % threshold of similarity
N = length(y); % length of time series
phi = zeros(2, 1); % phi(1)=phi_m, phi(2)=phi_{m+1}

for k = 1:2
	m = mnom + k - 1; % pattern length
	numVectors = N - (m - 1)*tau; % number of delay vectors of length m
	if numVectors < 2
		out = NaN; return % time series too short for this pattern length and delay
	end
	C = zeros(numVectors, 1);

	% Form delay vectors x from the time series y: x(i,:) = y(i:tau:i+(m-1)*tau).
	% Built via one vectorized indexing operation instead of a numVectors
	% iteration loop:
	idx = (1:numVectors)' + (0:m - 1)*tau;
	x = y(idx);

	for i = 1:numVectors
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
		C(i) = sum(dr) / numVectors; % Number of x(j) within r of x(i)
	end
	phi(k) = mean(log(C));
end
out = phi(1) - phi(2);

end
