function out = CO_RM_AMInformation(y,tau)
% CO_RM_AMInformation Automutual information (Rudy Moddemeijer implementation)
%
% Wrapper for Rudy Moddemeijer's information code to calculate automutual
% information.
%
%---INPUTS:
% y, the input time series
% tau, the time lag at which to calculate the automutual information

% ------------------------------------------------------------------------------
% Copyright (C) 2017, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

if nargin < 2 || isempty(tau)
    tau = 1; % Default is to calculate the automutual information at lag 1
end
if tau >= length(y)
    warning('Time series too short for a time lag of %u',tau);
    out = NaN; return
end

y1 = y(1:end-tau);
y2 = y(1+tau:end); % time-delayed version of y
out = RM_information(y1,y2);


end
