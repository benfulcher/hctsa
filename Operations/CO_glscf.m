function glscf = CO_glscf(y,alpha,beta,tau)
% CO_glscf      The generalized linear self-correlation function of a time series.
%
% This function was introduced in Queiros and Moyano in Physica A, Vol. 383, pp.
% 10--15 (2007) in the paper "Yet on statistical properties of traded volume:
% Correlation and mutual information at different value magnitudes"
%
% The function considers magnitude correlations.
%
%---INPUTS:
% y, the input time series
% alpha and beta are real and nonzero parameters
% tau is the time-delay (can also be 'tau' to set to first zero-crossing of the ACF)
%
% When alpha = beta estimates how values of the same order of magnitude are
% related in time
% When alpha ~= beta, estimates correlations between different magnitudes of the
% time series.

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

% ------------------------------------------------------------------------------
%% Check inputs and set defaults
% ------------------------------------------------------------------------------
if nargin < 4 || isempty(tau)
    tau = 'tau';
end

% Set tau to first zero-crossing of the autocorrelation function with the input 'tau'
if strcmp(tau,'tau')
    tau = CO_FirstZero(y,'ac');
end

y1 = abs(y(1:end-tau)); % take magnitudes
y2 = abs(y(1+tau:end)); % take magnitudes

glscf = (mean((y1.^alpha).*(y2.^beta)) - mean(y1.^alpha)*mean(y2.^beta)) / ...
     		    (sqrt(mean(y1.^(2*alpha)) - mean(y1.^alpha)^2) ...
     		          * sqrt(mean(y2.^(2*beta)) - mean(y2.^beta)^2));


end
