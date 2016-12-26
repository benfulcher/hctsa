function out = DN_Quantile(y,p)
% DN_Quantile   Quantile of the data vector
%
% Calculates the quantile value at a specified proportion, p, using the
% Statistics Toolbox function, quantile.
%
%---INPUTS:
% y, the input data vector
% p, the quantile proportion

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

if nargin < 2
    fprintf(1,'Using quantile p = 0.5 (median) by default\n');
    p = 0.5;
end
if ~isnumeric(p) || (p < 0) || (p > 1)
    error('p must specify a proportion, in (0,1)');
end

out = quantile(y,p);

end
