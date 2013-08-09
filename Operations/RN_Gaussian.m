% RN_Gaussian
% 
% Returns a random number drawn from a normal distribution.
% 
% INPUT:
% y, the input time series (but this input isn't used so it could be anything)
% 
% This operation was sometimes used as a control, as a kind of 'null operation'
% with which to compare to the performance of other operations that use
% informative properties of the time-series data, but is otherwise useless for
% any real analysis.
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

function out = RN_Gaussian(y)
% Ben Fulcher, 2009

if nargin < 1
    y = [];
end

out = randn(1,1);

end