% CO_FirstZero
% 
% Returns the first zero-crossing of a given autocorrelation function.
% 
% y, the input time series
% corrfn, the self-correlation function to measure:
%         (i) 'ac': normal linear autocorrelation function. Uses CO_AutoCorr to
%                   calculate autocorrelations.
% maxtau, a maximum time-delay to search up to
% 
% In future, could add an option to return the point at which the function
% crosses the axis, rather than the first integer lag at which it has already
% crossed (what is currently implemented)
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = CO_FirstZero(y,corrfn,maxtau)
% Ben Fulcher, 2008

N = length(y); % the length of the time series

if nargin < 2 || isempty(corrfn)
    corrfn = 'ac'; % autocorrelation by default
end
if nargin < 3 || isempty(maxtau)
    maxtau = N; % search up to a maximum of the length of the time series
    % maxtau = 400; % searches up to this maximum time lag
    % maxtau = min(maxtau,N); % searched up to the length of the time series if this is less than maxtau
end

% Select the self-correlation function as an inline function
% Eventually could add additional self-correlation functions
switch corrfn
case 'ac'
    Corr_fn = @(x) CO_AutoCorr(y,x); % autocorrelation at time lag x
otherwise
    error('Unknown correlation function ''%s''',corrfn);
end

% Calculate autocorrelation at increasing lags, until you find a negative one
for tau = 1:maxtau-1
    if Corr_fn(tau) < 0 % we know it starts positive (1), so first negative will be the zero-crossing
        out = tau; return
    end
end

% If haven't left yet, set output to maxtau
out = maxtau;

end