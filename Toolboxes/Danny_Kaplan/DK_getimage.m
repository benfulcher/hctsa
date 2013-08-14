% [data2, images] = DK_getimage(data,pred) finds the scalar images of
% the points in a time series <pred> time sets in the future
% data --- matrix of embedded data (from lagembed)
% pred --- look ahead time, default value 1
% Returns
% data2 --- a new embedded data matrix appropriately trimmed
% images --- the images (at time <pred>) of the points in data2
% This is a convenience program to trim an embedding appropriately.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 1996, D. Kaplan <kaplan@macalester.edu>
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

function [data2, images] = DK_getimage(data, pred)

if nargin < 2
  pred = 1;
end

images = data((1+pred):length(data),1);
data2 = data(1:(length(data)-pred),:);

end