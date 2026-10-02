function out = CO_FZCGLSCF(y, alpha, beta, maxtau)
% CO_FZCGLSCF   The first zero-crossing of the generalized self-correlation function.
%
% Returns the lag at which the generalized self-correlation function of Queiros and
% Moyano (see CO_GLSCF), the correlation between |y(t)|^alpha and |y(t+tau)|^beta,
% first changes sign as tau = 1, 2, ... increases. The crossing is placed by linear
% interpolation between the two lags either side. If the function never
% changes sign, the output is maxtau.
%
% ---INPUTS:
% y, the input time series
% alpha, the parameter alpha
% beta, the parameter beta
% maxtau, [optional] the maximum time delay to search up to (default: the length of
%         the time series)
%
% ---OUTPUTS:
% out, a scalar: the (interpolated) lag of the first zero-crossing, or maxtau.
%
% ---REFERENCES:
% Queiros and Moyano, "Yet on statistical properties of traded volume: Correlation and
% mutual information at different value magnitudes", Physica A 383, 10--15 (2007).

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

N = length(y); % the length of the time series

if nargin < 4 || isempty(maxtau)
	maxtau = N;
	% maxtau = 400; % searches up to this maximum time lag
	% maxtau = min(maxtau,N); % make sure no longer than the time series itself
end

glscfs = zeros(maxtau, 1);

for i = 1:maxtau
	tau = i;
	% y1 = abs(y(1:end-tau));
	% y2 = abs(y(1+tau:end));

	glscfs(i) = CO_GLSCF(y, alpha, beta, tau);

	if (i > 1) && (glscfs(i) * glscfs(i - 1) < 0)
		% Draw a straight line between these two and look at where hits zero
		out = i - 1 + glscfs(i - 1) / (glscfs(i - 1) - glscfs(i));
		return
	end
end

out = maxtau; % if the function hasn't exited yet, set output to maxtau

end
