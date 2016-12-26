function out = CO_tc3(y,tau)
% CO_tc3    Normalized nonlinear autocorrelation function, tc3.
%
% Computes the tc3 function, a normalized nonlinear autocorrelation, at a
% given time-delay, tau.
% Statistic is for two time-delays, normalized in terms of a single time delay.
% Used as a test statistic for higher order correlational moments in surrogate
% data analysis.
%
%---INPUTS:
% y, input time series
% tau, time lag
%
%---OUTPUTS:
% The raw tc3 expression, its magnitude, the numerator and its magnitude, and
% the denominator.
%
% See documentation of the TSTOOL package (http://www.physik3.gwdg.de/tstool/)
% for further details about this function
% (i.e., http://www.physik3.gwdg.de/tstool/manual.pdf)

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
% Set the time lag as a measure of the time-series correlation length
% ------------------------------------------------------------------------------
% Can set the time lag, tau, to be 'ac' or 'mi'
if strcmp(tau,'ac')
    tau = CO_FirstZero(y,'ac');
    % tau is first zero crossing of the autocorrelation function
elseif strcmp(tau,'mi')
    tau = CO_FirstMin(y,'mi');
    % tau is the first minimum of the automutual information function
end

% ------------------------------------------------------------------------------
% Compute tc3 statistic
% ------------------------------------------------------------------------------

yn = y(1:end-2*tau);
yn1 = y(1+tau:end-tau); % yn1, tau steps ahead
yn2 = y(1+2*tau:end); % yn2, 2*tau steps ahead

% The expression used in TSTOOL tc3:
out.raw = mean(yn.*yn1.*yn2)/abs(mean(yn.*yn1))^(3/2);

% The magnitude
out.abs = abs(out.raw);

% The numerator
out.num = mean(yn.*yn1.*yn2);
out.absnum = abs(out.num);

% The denominator
out.denom = abs(mean(yn.*yn1))^(3/2);

end
