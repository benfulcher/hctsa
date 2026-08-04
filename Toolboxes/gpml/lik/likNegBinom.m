function [varargout] = likNegBinom(link, hyp, y, mu, s2, inf, i)

% likNegBinom - Negative binomial likelihood function for count data y.
% The expression for the likelihood is 
%   likNegBinom(f) = 1/Z * mu^y / (r+mu)^(r+y), Z = G(y+1)*G(r)/(r^r*G(y+r))
% with G(t)=gamma(t)=(t-1)!, mean=mu and variance=mu*(mu+r)/r, where r is the 
% number of failures parameters, mu = g(f) is the negative binomial intensity,
% f is a Gaussian process and y is the non-negative integer count data.
% Hence, we have -- with log(Z)=r*log(r)+L(y+r)-L(y+1)-L(r), L(t)=gammaln(t) --
%   llik(f) = log(likNegBinom(f)) = y*log(mu)-(r+y)*log(r+mu)-log(Z).
%
% We provide two inverse link functions 'exp' and 'logistic':
%   g(f) = exp(f) and g(f) = log(1+exp(f))).
% The link functions are located at util/glm_invlink_*.m.
% 
% Note that for neither link function the likelihood lik(f) is log concave.
% Note further that for r->oo, we recover likPoisson.
% 
% The hyperparameters are:
%
% hyp = [  log(r)  ]
%
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likFunctions.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme.
%
% Copyright (c) by Hannes Nickisch, 2016-12-09.
%
% See also likFunctions.m.

if nargin<4, varargout = {'1'}; return; end   % report number of hyperparameters
if ~exist('psi'), mypsi = @digamma; else mypsi = @psi; end    % no psi in Octave

lr = hyp; r = exp(lr);

if nargin<6                              % prediction mode if inf is not present
  if numel(y)==0,  y = zeros(size(mu)); end
  s2zero = 1; if nargin>4&&numel(s2)>0&&norm(s2)>eps, s2zero = 0; end  % s2==0 ?
  if s2zero                                                    % log probability
    lZ = gammaln(y+1) + gammaln(r) - gammaln(y+r) - r*lr;
    lg = g(mu,link);
    mx = max(lg,lr); lgr = log(exp(lg-mx)+exp(lr-mx))+mx;       % log(exp(lg)+r)
    lp = y.*lg - (y+r).*lgr - lZ;
  else
    lp = likNegBinom(link, hyp, y, mu, s2, 'infEP');
  end
  ymu = {}; ys2 = {};
  if nargout>1                                 % compute y moments by quadrature
    n = max([length(y),length(mu),length(s2)]); on = ones(n,1);
    N = 20; [t,w] = gauher(N); oN = ones(1,N); lw = ones(n,1)*log(w');
    mu = mu(:).*on; sig = sqrt(s2(:)).*on;                        % vectors only
    lg = g(sig*t'+mu*oN,link); 
    ymu = exp(logsumexp2(lg+lw));     % first moment using Gaussian-Hermite quad
    if nargout>2
      elg = exp(lg);
      yv = elg.*(elg/r+1); % second y moment from negative binomial distribution
      ys2 = (yv+(elg-ymu*oN).^2)*w;
    end
  end
  varargout = {lp,ymu,ys2};
else
  switch inf 
  case 'infLaplace'
    [lg,dlg,d2lg,d3lg] = g(mu,link); elg = exp(lg);
    lZ = gammaln(y+1) + gammaln(r) - gammaln(y+r) - r*lr;
    mx = max(lg,lr); lgr = log(exp(lg-mx)+exp(lr-mx))+mx;       % log(exp(lg)+r)
    a = 1./(1+r./elg); da = a.*(1-a).*dlg;                % auxiliary quantities
    if nargin<7                                             % no derivative mode
      lp = y.*lg - (y+r).*lgr - lZ;
      dlp = {}; d2lp = {}; d3lp = {};                         % return arguments
      if nargout>1
        dlp = y.*dlg - (y+r).*a.*dlg;        % dlp, derivative of log likelihood
        if nargout>2                    % d2lp, 2nd derivative of log likelihood
          d2lp = y.*d2lg - (y+r).*(a.*d2lg + da.*dlg);
          if nargout>3                  % d3lp, 3rd derivative of log likelihood
            d3lp = y.*d3lg - (y+r).*(a.*d3lg + da.*(3*d2lg +(1-2*a).*dlg.*dlg));
          end
        end
      end
      varargout = {lp,dlp,d2lp,d3lp};
    else                                                       % derivative mode
      b = (y+r)./(elg+r);
      lp_dhyp = r*(1+log(r)-lgr-b-mypsi(r)+mypsi(y+r));
      dlp_dhyp = r*dlg.*a.*(b-1);                             % first derivative
      d2lp_dhyp = r*((d2lg.*a+dlg.*da).*(b-1)-(dlg.*a).^2.*b); % and also second
      varargout = {lp_dhyp,dlp_dhyp,d2lp_dhyp};
    end

  case 'infEP'
    if nargin<7                                             % no derivative mode
      % Since we are not aware of an analytical expression of the integral,
      % hence we use quadrature.
      varargout = cell(1,nargout);
      [varargout{:}] = lik_epquad({@likNegBinom,link},hyp,y,mu,s2);
    else                                                       % derivative mode
      varargout = {[]};                                     % deriv. wrt hyp.lik
    end

  case 'infVB'
    error('infVB not supported')
  end
end

% compute the log intensity using the inverse link function
function varargout = g(f,link)
  varargout = cell(nargout, 1);  % allocate the right number of output arguments
  if isequal(link,'exp')
    [varargout{:}] = glm_invlink_exp(f);
  elseif isequal(link,'logistic')
    [varargout{:}] = glm_invlink_logistic(f);
  else
    [varargout{:}] = glm_invlink_logistic2(link{2},f);
  end