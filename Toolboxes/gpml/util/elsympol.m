% Evaluate the order R elementary symmetric polynomial Newton's identity aka
% the Newton–Girard formulae: http://en.wikipedia.org/wiki/Newton's_identities
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-01-10.

function E = elsympol(Z,R)
% evaluate 'power sums' of the individual terms in Z
sz = size(Z); P = zeros([sz(1:2),R]); for r=1:R, P(:,:,r) = sum(Z.^r,3); end
E = zeros([sz(1:2),R+1]);                   % E(:,:,r+1) yields polynomial r
E(:,:,1) = ones(sz(1:2)); if R==0, return, end  % init recursion
E(:,:,2) = P(:,:,1);      if R==1, return, end  % init recursion
for r=2:R
  for i=1:r
    E(:,:,r+1) = E(:,:,r+1) + P(:,:,i).*E(:,:,r+1-i)*(-1)^(i-1)/r;
  end
end

