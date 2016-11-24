function out = CO_trev(y,tau)
% CO_trev   The trev function of a time series.
%
% Calculates the trev function, a normalized nonlinear autocorrelation,
% mentioned in the documentation of the TSTOOL nonlinear time-series analysis
% package (available here: http://www.physik3.gwdg.de/tstool/).
%
% The quantity is often used as a nonlinearity statistic in surrogate data
% analysis, cf. "Surrogate time series", T. Schreiber and A. Schmitz, Physica D,
% 142(3-4) 346 (2000).
%
%---INPUTS:
%
% y, time series
%
% tau, time lag (can be 'ac' or 'mi' to set as the first zero-crossing of the
%       autocorrelation function, or the first minimum of the automutual
%       information function, respectively)
%
%---OUTPUTS:
% the raw trev expression, its magnitude, the numerator and its magnitude, and
% the denominator.

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
%% Set defaults:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(tau)
    tau = 'ac';
end

% ------------------------------------------------------------------------------
% Can set the time lag, tau, to be 'ac' or 'mi'
% ------------------------------------------------------------------------------
if strcmp(tau,'ac')
    tau = CO_FirstZero(y,'ac');
    % tau is first zero crossing of the autocorrelation function
elseif strcmp(tau,'mi')
    tau = CO_FirstMin(y,'mi');
    % tau is the first minimum of the automutual information function
end

% ------------------------------------------------------------------------------
% Compute trev quantities
% ------------------------------------------------------------------------------

yn = y(1:end-tau);
yn1 = y(1+tau:end); % yn, tau steps ahead

% The trev expression used in TSTOOL:
out.raw = mean((yn1-yn).^3)/(mean((yn1-yn).^2))^(3/2);

% The magnitude
out.abs = abs(out.raw);

% The numerator
out.num = mean((yn1-yn).^3);
out.absnum = abs(out.num);

% The denominator
out.denom = (mean((yn1-yn).^2))^(3/2);

end
