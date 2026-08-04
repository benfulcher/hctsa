function [K,dK] = covSM(Q, hyp, x, z)

% Gaussian Spectral Mixture covariance function. The covariance function 
% parametrization depends on the sign of Q.
%
% Let t(Dx1) be an offset vector in dataspace e.g. t = x-z. Then w(DxP)
% are the weights and m(Dx|Q|) = 1/p, v(Dx|Q|) = (2*pi*ell)^-2 are spectral
% means (frequencies) and variances, where p is the period and ell the length
% scale of the Gabor function h(t2v,tm) given by the expression
%   h(t2v,tm) = exp(-2*pi^2*t2v).*cos(2*pi*tm)
%
% Then, the two covariances are obtained as follows:
%
% SM, spectral mixture:          Q>0 => P = 1
%   k(x,z) = w'*h((t.*t)'*v,t'*m), t = x-z
%
% SMP, spectral mixture product: Q<0 => P = D
%   k(x,z) = prod(w'*h(T*T*v,T*m)), T = diag(t), t = x-z
%
% Note that for D=1, the two modes +Q and -Q are exactly the same.
%
% The hyperparameters are:
%
% hyp = [ log(w(:))
%         log(m(:))
%         log(sqrt(v(:))) ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Note that the spectral density H(s) = F[ h(t) ] of covGaboriso is given by
% H(s) = N(s|m,v)/2 + N(s|-m,v)/2 where m=1/p is the mean and v=(2*pi*ell)^-2
% is the variance of a symmetric Gaussian mixture. Hence the covGaboriso
% covariance forms a basis for the class of stationary covariances since a
% weighted sum of covGaboriso covariances corresponds to an isotropic
% location-scale mixture of a symmetric Gaussian mixture in the spectral domain.
%
% For more details, see 
% [1] SM: Gaussian Process Kernels for Pattern Discovery and Extrapolation,
%     ICML, 2013, by Andrew Gordon Wilson and Ryan Prescott Adams,
% [2] SMP: GPatt: Fast Multidimensional Pattern Extrapolation with GPs,
%     arXiv 1310.5288, 2013, by Andrew Gordon Wilson, Elad Gilboa, 
%     Arye Nehorai and John P. Cunningham, and
% [3] Covariance kernels for fast automatic pattern discovery and extrapolation
%     with Gaussian processes, Andrew Gordon Wilson, PhD Thesis, January 2014.
%     http://www.cs.cmu.edu/~andrewgw/andrewgwthesis.pdf
% [4] http://www.cs.cmu.edu/~andrewgw/pattern/.
%
% For Q>0, covSM corresponds to Eq. 4.7 in Ref [3].
% For Q<0, covSM corresponds to Eq. 14 in Ref [2] and Eq. 5.3 in Ref [3] but our
% w here corresponds to  w^2 in Eq. 14.
%
% Copyright (c) by Andrew Gordon Wilson and Hannes Nickisch, 2016-05-06.
%
% See also covFunctions.m, covGaboriso.m, covGaborard.m.

if nargin<1, error('You need to provide Q.'), end
smp = Q<0; Q = abs(Q);                    % switch between covSM and covSMP mode
if nargin<3                                            % report no of parameters
  if smp, K = '3*D*'; else K = '(1+2*D)*'; end, K = [K,sprintf('%d',Q)]; return
end
if nargin<4, z = []; end                                   % make sure, z exists

% Note that we have two implementations covSMgabor and covSMfast. The former
% constructs a weighted sum of products of covGabor covariances using covMask,
% covProd, covScale and covSum while the latter is a standalone direct
% implementation. The latter tends to be faster.
if nargout > 1
  if smp
    [K,dK] = covSMgabor(Q,hyp,x,z,smp);
  else
    [K,dK] = covSMfast(Q,hyp,x,z,smp);               % faster direct alternative
  end
else
  K = covSMfast(Q,hyp,x,z,smp);                      % faster direct alternative
end

function [K,dK] = covSMgabor(Q,hyp,x,z,smp)
  [n,D] = size(x); P = smp*D+(1-smp);               % dimensionality, P=D or P=1
  lw = reshape(hyp(         1:P*Q) ,P,Q);                  % log mixture weights
  lm = reshape(hyp(P*Q+    (1:D*Q)),D,Q);                   % log spectral means
  ls = reshape(hyp(P*Q+D*Q+(1:D*Q)),D,Q);     % log spectral standard deviations
  if smp % 1) the product of weighted sums of 1d covGabor functions or
    fac = cell(1,D);
    for d=1:D
      add = cell(1,Q);       % a) addends for weighted sum of 1d Gabor functions
      for q=1:Q, add{q} = {'covScale',{'covMask',{d,{'covGaboriso'}}}}; end
      fac{d} = {'covSum',add};                     % b) combine addends into sum
    end
    fac{D+1} = 0;                       % disable cache to avoid memory problems
    cov = {'covProd',fac};                     % c) combine factors into product
  else   % 2) the weighted sum of multivariate covGaborard covariance functions.
                                  % weighted sum of multivariate Gabor functions
    add = cell(1,Q); for q=1:Q, add{q} = {'covScale',{'covGaborard'}};  end
    cov = {'covSum',add};                                     % combine into sum
  end
  if smp    % assemble hyp; covGabor is parametrised using -ls-log(2*pi) and -lm
    hypgb = [-log(2*pi)-ls(:)'; -lm(:)'; lw(:)'/2];
  else
    hypgb = [-log(2*pi)-ls;     -lm;     lw/2    ];
  end
  if nargout>1
    [K,dK] = feval(cov{:},hypgb(:),x,z);
    dK = @(R) dirder_gabor(R,dK,Q,x,smp);
  else
    K = feval(cov{:},hypgb(:),x,z);
  end
  
function [dhyp,dx] = dirder_gabor(R,dK,Q,x,smp)
  if nargout > 1
    [dhyp,dx] = dK(R);
  else
    dhyp = dK(R);
  end
  D = size(x,2);
  if smp
    dhyp = reshape(dhyp,3,D*Q);
    dls = -dhyp(1,:); dlm = -dhyp(2,:); dlw = 0.5*dhyp(3,:);
  else
    dhyp = reshape(dhyp,2*D+1,Q);
    dls = -dhyp(1:D,:); dlm = -dhyp(D+1:2*D,:); dlw = 0.5*dhyp(2*D+1,:);
  end
  dhyp = [dlw(:); dlm(:); dls(:)];

function [K,dK] = covSMfast(Q,hyp,x,z,smp)
  xeqz = isempty(z); dg = strcmp(z,'diag');           % sort out different types
  [n,D] = size(x); P = smp*D+(1-smp);               % dimensionality, P=D or P=1
  w = exp(reshape(  hyp(         1:P*Q) ,P,Q));                % mixture weights
  m = exp(reshape(  hyp(P*Q+    (1:D*Q)),D,Q));                 % spectral means
  v = exp(reshape(2*hyp(P*Q+D*Q+(1:D*Q)),D,Q));             % spectral variances
  if dg
    T = zeros(n,1,D);
  else
    if xeqz
      T = 2*pi*bsxfun(@minus,reshape(x,n,1,D),reshape(x,1,n,D));
    else
      T = 2*pi*bsxfun(@minus,reshape(x,n,1,D),reshape(z,1,[],D));
    end
  end, T = reshape(T,[],D);
  if smp
    h = @(t2v,tm) exp(-0.5*t2v).*cos(tm);                       % Gabor function
    K = 1; w = reshape(w,Q,P)'; m = reshape(m,Q,D)'; v = reshape(v,Q,D)';
    for d=1:D
      K = K .* ( h( (T(:,d).*T(:,d))*v(d,:), T(:,d)*m(d,:) )*w(d,:)' );
    end
    K = reshape(K.*ones(size(T,1),1),n,[]);
  else
    E = exp(-0.5*(T.*T)*v); H = E.*cos(T*m);
    K = reshape(H*w',n,[]);
    if nargout>1
      vec = @(x) x(:);   
      dKdhyp = @(R) [ (H'*R(:)).*w'; 
                     -vec( ((E.*sin(T*m).*(R(:)*w))'* T    )'.*m );
                     -vec( ((H.*          (R(:)*w))'*(T.*T))'.*v ) ];
      dK = @(R) dirder_fast(R,E,T,H,dKdhyp,x,z, m,v,w);
    end
  end
  
function [dhyp,dx] = dirder_fast(R, E,T,H,dKdhyp,x,z, m,v,w)
  dhyp = dKdhyp(R);
  if nargout>1
    xeqz = isempty(z); dg = strcmp(z,'diag'); [n,D] = size(x);
    if dg
      dx = zeros(size(x));
    else
      A = reshape( (R(:)*w).*E.*sin(T*m)*m' + (((R(:)*w).*H)*v').*T, n, [], D);
      dx = -2*pi*squeeze(sum(A,2));
      if xeqz, dx = dx + 2*pi*squeeze(sum(A,1)); end
    end
  end