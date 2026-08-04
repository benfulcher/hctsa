function [K,dK] = covW(i, hyp, x, z)
% covW i-times integrated Wiener process covariance function.
%
%   [K,dK] = covW(i, hyp, x, z)
%
% For i=-1, this is just the white noise covariance, see covNoise.
% For i= 0, this is the Wiener process covariance,
% for i= 1, this is the integrated Wiener process covariance (velocity),
% for i= 2, this is the twice-integrated Wiener process covariance (accel.),
% for i= 3, this is the thrice-integrated Wiener process covariance.
%
% If dw(x) is a Wiener process then dw(x+s^2)-dw(x) = N(0,s^2).
% The Wiener process is the integral of a white noise process, see covEye.
%
% The hyperparameters are:
%
%  hyp = [ log(sf) ]
%
% where SF is the amplitude of the covariance, see below.
%
% The covariance function -- given that x,z>=0 -- is specified as:
%
% k(x,z) = sf^2 * ki(x,z), where ki(x,z) is given for
% different values of i by the following expressions
% i = -1, ki(x,z) =  \delta(x,z)
% i >= 0, ki(x,z) = 1/ai*min(x,z)^(2*i+1) + bi*min(x,z)^(i+1) * |x-z| * ri(x,z),
% with the coefficients ai, bi and the residual ri(x,z) defined as follows:
% i = 0, ai =   1, bi = 0
% i = 1, ai =   3, bi = 1/  2, ri(x,z) = 1
% i = 2, ai =  20, bi = 1/ 12, ri(x,z) = x+z-1/2*min(x,z)
% i = 3, ai = 252, bi = 1/720, ri(x,z) = 5*max(x,z)^2+2*x*z+3*min(x,z)^2
%
% For more help on design of covariance functions, try "help covFunctions".
%
% References:
% See the paper Probabilistic ODE Solvers with Runge-Kutta Means by Schober,
% Duvenaud and Hennig, NIPS, 2014, for more details.
%
% See also cov/covEye.m, covFunctions.m.
%
% Copyright (c) by Hannes Nickisch, 2017-10-05.

if nargin<3, K = '1'; return; end                  % report number of parameters
if nargin<4, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

[n,D] = size(x); ox = ones(n,1);
if D~=1, error('Covariance is defined for 1d data only.'), end
if any(x<0) || (~xeqz && any(z<0))
  error('Covariance is defined for nonnegative data only.')
end
sf = exp(hyp(1));

if dg                                                               % vector kxx
  Mn = x; Mx = x; I = ox; D = zeros(n,1); S = 2*x; P = x.*x;
else
  if xeqz                                                 % symmetric matrix Kxx
    I = bsxfun(@le,x,x'); Mn = ox*x'; T = x*ox'; P = x*x';
  else                                                   % cross covariances Kxz
    oz = ones(size(z)); I = bsxfun(@le,x,z'); Mn = ox*z'; T = x*oz'; P = x*z';
  end
  D = T-Mn; S = T+Mn;                      % D = x-z, S = x+z, P = x*z, I = x<=z
  Mx = Mn; Mx(~I) = T(~I); Mn(I) = T(I);          % Mn = min(x,z), Mx = max(x,z)
end

switch i
  case -1, K = double(abs(D)<eps*eps);
  case  0, K = Mn;
  case  1, K = Mn.^3/  3 + abs(D).*Mn.^2/  2;
  case  2, K = Mn.^5/ 20 + abs(D).*Mn.^3/ 12.*(S-Mn/2);
  case  3, K = Mn.^7/252 + abs(D).*Mn.^4/720.*(5*Mx.^2+2*P+3*Mn.^2);
  otherwise, error('unknown degree')
end
K = sf^2*K;                                                    % signal variance

if nargout > 1
  dK = @(Q) dirder(Q,K,I,Mn,Mx,D,S,P,sf,i,dg,xeqz,x);        % directional deriv
end

function [dhyp,dx] = dirder(Q,K,I,Mn,Mx,D,S,P,sf,i,dg,xeqz,x)
  dhyp = 2*(Q(:)'*K(:));                                       % signal variance
  switch i
    case -1, dMn = 0; dD = 0; dS = 0;
    case  0, dMn = 1; dD = 0; dS = 0;
    case  1, dMn = Mn.^2+abs(D).*Mn; dD = Q.*sign(D).*Mn.^2/2; dS = 0;
    case  2, dMn = Mn.^4/4 + abs(D).*Mn.^2.*(S/4-Mn/6);
             dD = Q.*sign(D).*Mn.^3/12.*(S-Mn/2); dS = Q.*abs(D).*Mn.^3/12;
    case  3, dMn = Mn.^6/36 + abs(D).*Mn.^3/180.*(5*Mx.^2+2*P+9/2*Mn.^2);
             dD = Q.*sign(D).*Mn.^4/720.*(5*Mx.^2+2*P+3*Mn.^2);
             dS = 0;                       % dMx and dP are assumend to be small
  end
  Q = Q.*dMn;
  if nargout > 1
    if dg
      dx = zeros(size(x));
    else
      if xeqz
        dx = sum(Q.*I,2)+sum(Q.*(1-I),1)'+sum(dD,2)-sum(dD,1)'+sum(dS+dS',2);
      else
        dx = sum(Q.*I,2)+sum(dD,2)+sum(dS,2);
      end
    end
    dx = sf^2*dx;                                              % signal variance
  end
