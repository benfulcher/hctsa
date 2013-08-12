function d = DK_disttyp(z,percs)
% d = disttyp(z,percs)
% Calculates typical distances between pre-images
% z - embedded data
% percs -- percentiles to use
% returns distances at the given percentiles
% all of this is from a small sample of all pairs of distances
% Copyright (c) 1996 by D. Kaplan, All Rights Reserved

if any(percs< 0 | percs > 1)
    error('disttyp: percs must be in [0,1]');
  end
Npts = length(z);
% we want up to 10000 points
ntargs = ceil(10000/(Npts-1));
ntargs = min(Npts,ntargs);
dists = zeros(ntargs*Npts,1);
count=0;

for j=1:ntargs
  dists((count+1):(count+Npts)) = DK_onedist(z,z(j,:));
  count = count+Npts;
end

dists = sort(dists(1:count));
% get rid of the zeros corresponding to self-distances
dists = dists((ntargs+1):length(dists) );

d = zeros(size(percs));

for p = 1:(length(percs))
  ind = round(percs(p)*(length(dists)-1));
  ind = max(ind,0);
  ind = min(ind,length(dists)-1);
  d(p) = dists( ind+1 );
end

