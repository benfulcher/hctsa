% DK_findneib(z, pt, k, r ) finds the nearest neighbors to pt in z
% z -- matrix of points, 1 per row
% pt -- vector of a single point
% k -- number to find
% r (optional -- if specified, find all neighbors closer than this
%
% inds -- indices of the closest points to pt
% dist -- corresponding distances from pt
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

function [inds, dist] = DK_findneib(z, pt, k, r )

ds = DK_onedist(z,pt);
[Y,I] = sort(ds);
if nargin == 3 
  inds = I( 1:k );
  dist = Y( 1:k );
else
  foo = find(Y<=r);
  inds = I(foo);
  dist = Y(foo);
end


