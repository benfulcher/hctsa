% CO_AutoCorr
% 
% Computes the autocorrelation of an input time series, y, at a time-lag, tau
% 
% INPUTS:
% y, a scalar time series column vector
% tau, the time-delay. If tau is a scalar, returns autocorrelation for y at that
%       lag. If tau is a vector, returns autocorrelations for y at that set of
%       lags.
%       
% Output is the autocorrelation at the given time-lag
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

function out = CO_AutoCorr(y,tau)

% Check inputs:
if nargin < 2 || isempty(tau)
    tau = 1;
end
N = length(y); % length of the time sries

if length(tau) == 1 % output a single value at the given time-lag
    out = sum((y(1:N-tau) - mean(y(1:N-tau))).*(y(tau+1:N) ...
	            - mean(y(tau+1:N))))/N/std(y(1:N-tau))/std(y(tau+1:N));
else % output values over a range of time-lags
    out = zeros(length(tau),1);
    for i = 1:length(tau)
        out(i) = sum((y(1:N-tau(i)) - mean(y(1:N-tau(i)))).*(y(tau(i)+1:N) ...
			        - mean(y(tau(i)+1:N))))/N/std(y(1:N-tau(i)))/std(y(tau(i)+1:N));
    end
end

end