function out = CO_fzcglscf(y,alpha,beta,maxtau)
% CO_fzcglscf   The first zero-croessing of the generalized self-correlation function
%
% Returns the first zero-crossing of the generalized self-correlation function
% introduced in Duarte Queiros and Moyano in Physica A, Vol. 383, pp. 10--15
% (2007) in the paper "Yet on statistical properties of traded volume:
% Correlation and mutual information at different value magnitudes"
% Uses CO_glscf to calculate the generalized self-correlations.
% Keeps calculating until the function finds a minimum, and returns this lag.
%
%---INPUTS:
% y, the input time series.
% alpha, the parameter alpha.
% beta, the parameter beta.
% maxtau [opt], a maximum time delay to search up to (default is the time-series
%                length).

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
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

glscfs = zeros(maxtau,1);

for i = 1:maxtau
	tau = i;
    % y1 = abs(y(1:end-tau));
    % y2 = abs(y(1+tau:end));

    glscfs(i) = CO_glscf(y,alpha,beta,tau);

	if (i > 1) && (glscfs(i)*glscfs(i-1) < 0)
		% Draw a straight line between these two and look at where hits zero
		out = i - 1 + glscfs(i)/(glscfs(i)-glscfs(i-1));
		return
	end
end

out = maxtau; % if the function hasn't exited yet, set output to maxtau

end
