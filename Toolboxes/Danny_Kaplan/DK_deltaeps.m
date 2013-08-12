function [delta,epsilon] = DK_deltaeps(z, images, lockout)
% [delta,epsilon] = deltaeps(z, images)
% Delta-epsilon method
% z      -- embedded data as from getimage()
% images -- as from getimage()
% lockout-- don't consider points closer in time than this
% --------------
% delta  -- distances between pre-images
% eps    -- distances between corresponding images
% plot(delta,eps,'.') shows how image distance depends on 
% pre-image distance
% Copyright (c) 1996 by D. Kaplan, All Rights Reserved

if( nargin < 3 )
  lockout = 0;
end
npts = length(z);
k=min(npts-2,20);
delta=zeros(k*npts,1);
epsilon=delta;
count=0;

for i = 1:npts
  % find the nearest neighbors
  [inds,dists] = DK_findneib(z, z(i,:), k+10);
  % eliminate the ones closer in time than the lockout
  foo = find( abs(inds - i) > lockout );
  inds = inds(foo);
  dists = dists(foo);
  kk = min(length(dists),k);
  % store the distances between pre-images in the output
  delta((count+1):(count+kk),:) = dists(1:kk);
  % store the distances between images
  epsilon((count+1):(count+kk),:) = images(inds(1:kk)) - images(i);
  count = count+kk;
end
 
% take the absolute value of the distances between images
epsilon = abs(epsilon(1:count));  
delta = delta(1:count);