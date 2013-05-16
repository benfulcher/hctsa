function res = KP_onedist(z,pt)
% ONEDIST(z,pt) calculates the distance between point pt and each
% row in matrix z
% Copyright (c) 1996 by D. Kaplan, All Rights Reserved

[r,c] = size(z);
if c ~= length(pt)
  error('pt and z must have same number of columns');
end

sum = zeros(r,1);
for n=1:c
  foo = z(:,n) - pt(n);
  sum = sum + foo.*foo;
end
res = sqrt(sum);
