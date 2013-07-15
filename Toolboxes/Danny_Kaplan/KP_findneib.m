function [inds, dist] = KP_findneib(z, pt, k, r )
% FINDNEIB(z, pt, k, r ) finds the nearest neighbors to pt in z
% z -- matrix of points, 1 per row
% pt -- vector of a single point
% k -- number to find
% r (optional -- if specified, find all neighbors closer than this
%
% inds -- indices of the closest points to pt
% dist -- corresponding distances from pt
% Copyright (c) 1996 by D. Kaplan, All Rights Reserved

ds = KP_onedist(z,pt);
[Y,I] = sort(ds);
if nargin == 3 
  inds = I( 1:k );
  dist = Y( 1:k );
else
  foo = find(Y<=r);
  inds = I(foo);
  dist = Y(foo);
end


