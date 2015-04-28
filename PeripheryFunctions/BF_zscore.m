% ------------------------------------------------------------------------------
% BF_zscore
% ------------------------------------------------------------------------------
% 
% Applies a z-score to the input without using a Statistics Toolbox licence.
% 
%---INPUT:
% x, the input time series (or any vector).
% 
%---OUTPUT:
% z, the z-scored transformation of the input.
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

function zscoredData = BF_zscore(inputData)

% By default, z-score twice to reduce the numerical error:
zscoredData = (inputData - mean(inputData)) / std(inputData);
zscoredData = (zscoredData - mean(zscoredData)) / std(zscoredData);

end