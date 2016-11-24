function colormap = BF_MakeBrightenedColorMap(rgb,numGrads)
% BF_MakeBrightenedColorMap  Given an rgb vector, progressively brightens into a colormap
%
%---INPUTS:
% xy, a vector (or cell) of x-y co-ordinates of points on the plot
% TimeSeries, a structure array of time series making up the plot
% annotateParams, structure of custom plotting parameters

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
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
    numGrads = 8;
end

colormap = zeros(numGrads,3);

brightenRange = linspace(0,1,numGrads);
for i = 1:numGrads
    colormap(i,:) = brighten(rgb,brightenRange(i));
end

end
