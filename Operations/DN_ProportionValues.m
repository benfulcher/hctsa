function out = DN_ProportionValues(x,propWhat)
% DN_ProportionValues   Proportion of values in a data vector.
%
% Returns statistics on the values of the data vector: the proportion of zeros,
% the proportion of positive values, and the proportion of values greater than or
% equal to zero.
%
%---INPUTS:
% x, the input time series
%
% propWhat, the proportion of a given type of value in the time series:
%           (i) 'zeros': values that equal zero
%           (ii) 'positive': values that are strictly positive
%           (iii) 'geq0': values that are greater than or equal to zero

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

%-------------------------------------------------------------------------------
% Check Inputs:
%-------------------------------------------------------------------------------
if nargin < 2
    propWhat = 'positive';
end

N = length(x); % length of the time series

switch propWhat
    case 'zeros' % returns the proportion of zeros in the input vector
        out = sum(x == 0)/N;

    case 'positive'
        out = sum(x > 0)/N;

    case 'geq0'
        out = sum(x >= 0)/N;

    otherwise
        error('Unknown condition to measure: ''%s''',propWhat);
end

end
