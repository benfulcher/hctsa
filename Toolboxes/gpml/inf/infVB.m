function [post, nlZ, dnlZ] = infVB(hyp, mean, cov, lik, x, y)

% Variational approximation to the posterior Gaussian process with MKL 
% covariance function hyperparameter optimisation.
% The function takes a likelihood function (see likFunction.m), and is designed
% to be used with gp.m. See also infFunctions.m.
%
% Minimisation of an upper bound on the negative marginal likelihood 
% \int N(f|0,K) p(y|f) df \ge Psi(ga,theta) with hyperparameters 
% theta=cov.hyp and likelihood lik(y,f) = p(y|f).
%
% Psi(ga,theta) = nlZ =
%     = (ln|K+ga|-ln|ga| + h(ga) - b'*inv(inv(K)+inv(ga)))*b)/2
%     = (ln|A| + ln|K|   + h(ga) - b'*inv(A)*b)/2 with A = inv(K)+inv(ga)
%
% The map (ga,K) |-> Psi(ga,theta) - ln|K| is jointly convex and ln|K| is
% concave in its linear parameters theta.
%
% We optimise the convex Psi(ga,theta)-ln|K|+z'*theta/2  +  B(t,theta) criterion
% in the outer loop using the barrier function B(t,theta) = -sum(log(theta))/t 
% to enforce positivity w.r.t. theta. 
% The linear term z'*theta/2 is an upper bound on the concave ln|K| term.
% We use an interleaved optimisation with inner loops doing Newton steps in ga
% and an outer loop doing joint Newton updates in (ga,theta). Line searches are
% done using derivative-free 'Brent's minimum' search.
%
% Copyright (c) by Hannes Nickisch 2012-11-07.
%
% See also INFMETHODS.M.

maxitinner = 20;                            % number of inner Newton steps in ga

% some less important parameters
ep = 1e-8;    % small constant for the ridge added to the stab. Newton direction
smax = 5; Nline = 15; thr = 1e-4;                       % line search parameters

tol = 1e-7;                   % tolerance for when to stop the Newton iterations
inf = 'infVB';
[n,D] = size(x);
K = feval(cov{:}, hyp.cov, x);                  % evaluate the covariance matrix
m = feval(mean{:}, hyp.mean, x);                      % evaluate the mean vector
if iscell(lik), likstr = lik{1}; else likstr = lik; end
if ~ischar(likstr), likstr = func2str(likstr); end

if    (norm(m)>1e-10 || numel(hyp.mean)>0) ...
   && (strcmp(likstr,'likErf')||strcmp(likstr,'likLogistic')) 
    error('only meanZero implemented for classification')
end
y = y-m;         % no we have either zero mean or non-classification likelihoods
if strcmp(likstr,'likGauss')
  ga = exp(2*hyp.lik)*ones(n,1);      % best ga is known for Gaussian likelihood
elseif strcmp(likstr,'likErf')
  ga = ones(n,1);        % use a fixed ga for errf likelihood due to asymptotics
else
  % INNER compute the Newton direction of ga
  ga = ones(n,1);                                               % initial values
  itinner = 0;
  nlZ_new = 1e100; nlZ_old = Inf;                  % make sure while loop starts
  while nlZ_old-nlZ_new>tol && itinner<maxitinner                 % begin Newton
    itinner = itinner+1;
    [nlZ_old,dga,d2ga] = Psi(ga,K,inf,hyp,lik,y);       % calculate grad/Hessian   
    ddga = -(d2ga+ep*eye(n))\dga;    % stabilized Newton direction + line search
    s = .99*min(ga./max(-ddga,0)); s = min(s,smax);    % max s, s.t. ga+s*ddga>0
    Psi_s = @(s,ddga,ga,K,inf,hyp,lik,y) Psi(ga+s*ddga,K,inf,hyp,lik,y);
    [s,nlZ_new] = brentmin(0,s,Nline,thr,   Psi_s,0,ddga,ga,K,inf,hyp,lik,y);
    ga = abs(ga + s*ddga);                       % update variational parameters
  end
end

[nlZ,dnlZ,d2nlZ,b] = Psi(ga,K,inf,hyp,lik,y);       % upp bd on neg log marg lik
W = 1./ga; sW = sqrt(W);                       % return the posterior parameters
L  = chol(eye(n)+sW*sW'.*K);                          % L'*L = B =eye(n)+sW*K*sW
iKtil = repmat(sW,1,n).*solve_chol(L,diag(sW));       % sW*B^-1*sW=inv(K+inv(W))
alpha = b - iKtil*(K*b);
post.alpha = alpha; post.sW = sW; post.L  = L;

if nargout>2                                           % do we want derivatives?
  dnlZ = hyp;                                   % allocate space for derivatives
  for j=1:length(hyp.cov)                                    % covariance hypers
    dK = feval(cov{:}, hyp.cov, x, [], j);
    if j==1, v = iKtil*(b./W); end
    dnlZ.cov(j) = sum(sum(iKtil.*dK))/2 - (v'*dK*v)/2; % implicit derivative = 0
  end
  if ~strcmp(likstr,'likGauss')                              % likelihood hypers
    for j=1:length(hyp.lik)
      dhhyp = feval(lik{:},hyp.lik,y,[],ga,inf,j);
      dnlZ.lik(j) = sum(dhhyp)/2;                      % implicit derivative = 0
    end
  else                                 % special treatment for the Gaussian case
    dnlZ.lik = sum(sum( (L'\eye(n)).^2 )) - exp(2*hyp.lik)*(alpha'*alpha);
  end
  for j=1:length(hyp.mean)                                         % mean hypers
    dm = feval(mean{:}, hyp.mean, x, j);
    dnlZ.mean(j) = -alpha'*dm;                         % implicit derivative = 0
  end
end

% variational lower bound along with derivatives w.r.t. ga and theta=hyp.cov
% psi = (ln|K+ga|-ln|ga| + h(ga) - b'*inv(inv(K)+inv(ga)))*b)/2
%     = (ln|A| + ln|K|   + h(ga) - b'*inv(A)*b)/2 with A = inv(K)+inv(ga)
% the code is numerically stable
function [nlZ,dga,d2ga,b] = Psi(ga,K,inf,hyp,lik,y)
  n = size(K,1);
  [h,b,dh,db,d2h,d2b] = feval(lik{:},hyp.lik,y,[],ga,inf);
  W = 1./ga; sW = sqrt(W);
  L = chol(eye(n)+sW*sW'.*K);                     % sum(log(diag(L))) = log|B|/2
  C = L'\(repmat(sW,1,n).*K);
  t = C*b;                         % t'*t-b'*K*b = -b'*inv(inv(K)+diag(1./ga))*b
  nlZ = sum(log(diag(L)))   +   ( sum(h)   +   t'*t-b'*K*b )/2;
  if nargout>1                     % Hessian w.r.t. variational parameters gamma
    iKtil = repmat(sW,1,n).*solve_chol(L,diag(sW)); % sW*inv(B)*sW=inv(K+inv(W))
    Khat = K-C'*C; v = Khat*b;                % K-K*sW*inv(B)*sW*K=inv(inv(K)+W)
    dga = ( diag(iKtil)-1./ga   +   dh   -   (v./ga).^2 - 2*v.*db )/2;
    if nargout>2                  % gradient w.r.t. variational parameters gamma
      w = v./ga.^2;
      d2ga = ( -iKtil.^2+diag(1./ga.^2)   +   diag(d2h) )/2 ...
              -Khat.*(w*(w+2*db)') -Khat.*(db*db') +diag(w.^2.*ga)-diag(v.*d2b);
      d2ga = (d2ga+d2ga')/2;                                        % symmetrise
    end
  end
