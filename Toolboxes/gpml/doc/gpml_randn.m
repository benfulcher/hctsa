function x = gpml_randn(seed,varargin)

% Generate pseudo-random numbers in a quick and dirty way.
% The function makes sure, we obtain the same random numbers using Octave and
% Matlab for the demo scripts.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2010-07-07.

if nargin==2
  sz = varargin{1};
else
  sz = zeros(1,nargin-1);
  for i=1:nargin-1
    sz(i) = varargin{i};
  end
end

n = prod(sz); N = ceil(n/2)*2;

% minimal uniform random number generator for uniform deviates from [0,1]
% by Park and Miller
a = 7^5; m = 2^31-1; 
% using Schrage's algorithm
q = fix(m/a); r = mod(m,a); % m = a*q+r
u = zeros(N+1,1); u(1) = fix(seed*2^31);
for i=2:N+1
  % Schrage's algorithm for mod(a*u(i),m)
  u(i) = a*mod(u(i-1),q) - r*fix(u(i-1)/q);
  if u(i)<0, u(i) = u(i)+m; end
end
u = u(2:N+1)/2^31;

% Box-Muller transform: Numerical Recipies, 2nd Edition, $7.2.8
% http://en.wikipedia.org/wiki/Box-Muller_transform                
w = sqrt(- 2*log(u(1:N/2)));                             % split into two groups
x = [w.*cos(2*pi*u(N/2+1:N)); w.*sin(2*pi*u(N/2+1:N))]; 
x = reshape(x(1:n),sz);
