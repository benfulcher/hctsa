% EN_Shannonpdf
% 
% Estimates the (log_2) Shannon entropy from the probability distribution of the time
% series.
% 
% INPUT:
% y, the input time series
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

function out = EN_Shannonpdf(y)
% Ben Fulcher, 2009

N = length(y); % time-series length
tmp = sort(y);
tmp = diff(tmp); % incremental differences of sorted set of values in y
if size(tmp,1) < size(tmp,2)
    tmp = tmp';
end

pdf = diff([0; find(tmp > 0); N])/N; % Empirical probability distribution
out = -sum(pdf.*log2(pdf));

end