% Evaluate the order R elementary symmetric polynomials using Newton's identity,
% the Newton-Girard formulae: http://en.wikipedia.org/wiki/Newton's_identities
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-01-10.
%                  speedup contributed by Truong X. Nghiem,   2016-01-20.

function E = elsympol(Z,R)
sz = size(Z);               % evaluate 'power sums' of the individual terms in Z
E = zeros([sz(1:2),R+1]);                       % E(:,:,r+1) yields polynomial r
% fast and efficient version of: for r=1:R, P(:,:,r) = sum(Z.^r,3); end
Zr = Z; P(:,:,1) = sum(Zr,3); for r=2:R, Zr = Zr.*Z; P(:,:,r) = sum(Zr,3); end
E(:,:,1) = ones(sz(1:2)); if R==0, return, end                  % init recursion
E(:,:,2) = P(:,:,1);      if R==1, return, end                  % init recursion
for r=2:R
  E(:,:,r+1) = P(:,:,1).*E(:,:,r)/r;              % i=1 is simpler than the rest
  for i=2:r
    E(:,:,r+1) = E(:,:,r+1) + P(:,:,i).*E(:,:,r+1-i)*(-1)^(i-1)/r;
  end
end