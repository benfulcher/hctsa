function out = CO_FirstZero(y,corrFun,maxTau)
% CO_FirstZero      The first zero-crossing of a given autocorrelation function
%
%---INPUTS:
%
% y, the input time series
% corrFun, the self-correlation function to measure:
%         (i) 'ac': normal linear autocorrelation function. Uses CO_AutoCorr to
%                   calculate autocorrelations.
% maxTau, a maximum time-delay to search up to.
%
% In future, could add an option to return the point at which the function
% crosses the axis, rather than the first integer lag at which it has already
% crossed (what is currently implemented).

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
% Check inputs:
% ------------------------------------------------------------------------------

N = length(y); % the length of the time series

if nargin < 2 || isempty(corrFun)
    corrFun = 'ac'; % autocorrelation by default
end
if nargin < 3 || isempty(maxTau)
    maxTau = N; % search up to a maximum of the length of the time series
    % maxTau = 400; % searches up to this maximum time lag
    % maxTau = min(maxTau,N); % searched up to the length of the time series if this is less than maxTau
end

% ------------------------------------------------------------------------------
% Select the self-correlation function as an inline function
% Eventually could add additional self-correlation functions
switch corrFun
case 'ac'
    % Autocorrelation at all time lags
    corrs = CO_AutoCorr(y,[],'Fourier');
    corrs = corrs(2:end); % remove the zero-lag result
otherwise
    error('Unknown correlation function ''%s''',corrFun);
end

% Calculate autocorrelation at increasing lags, until you find a negative one
for tau = 1:maxTau-1
    if corrs(tau) < 0 % we know it starts positive (1), so first negative will be the zero-crossing
        out = tau; return
    end
end

% If haven't left yet, set output to maxTau
out = maxTau;

end
