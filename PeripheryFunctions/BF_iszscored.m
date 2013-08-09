% BF_iszscored
% 
% Checks whether the input time series is (eps-close to being) z-scored.
% Used for displaying warning messages for functions that require z-scored
% inputs.
% 
% INPUT:
% x, the input time series (or any vector)
% 
% OUTPUT:
% iszscored, a logical with the verdict.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
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

function iszscored = BF_iszscored(x)
% Ben Fulcher, 2013

iszscored = ((abs(mean(x)) < eps) || (abs(std(x)-1) < eps));

end