function out = CO_RM_AMInformation(y,tau)
% CO_RM_AMInformation
%
% Wrapper for Rudy Moddemeijer's information code to calculate automutual
% information.
%
%---INPUTS:
% y, the input time series
% tau, the time lag at which to calculate the automutual information

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

if nargin < 2 || isempty(tau)
    tau = 1; % Default is to calculate the automutual information at lag 1
end

y1 = y(1:end-tau);
y2 = y(1+tau:end); % time-delayed version of y
out = RM_information(y1,y2);


end
