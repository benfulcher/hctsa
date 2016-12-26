function out = DN_Burstiness(y)
% DN_Burstiness     Burstiness statistic of a time series
%
% Returns the 'burstiness' statistic from
%
% Goh and Barabasi, 'Burstiness and memory in complex systems' Europhys. Lett.
% 81, 48002 (2008).
%
%---INPUT:
% y, the input time series
%
%---OUTPUT:
% The burstiness statistic, B.

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

r = std(y)/mean(y); % coefficient of variation

%-------------------------------------------------------------------------------
% Original Goh and Barabasi burstiness statistic, B:
out.B = (r - 1)/(r + 1);
% B = (std(y) - mean(y))/(std(y) + mean(y));

%-------------------------------------------------------------------------------
% Improved burstiness statistic, accounting for scaling for finite time series
% Kim and Jo, 2016, http://arxiv.org/pdf/1604.01125v1.pdf
N = length(y);
out.B_Kim = (sqrt(N+1)*r - sqrt(N-1))/((sqrt(N+1)-2)*r + sqrt(N-1));

end
