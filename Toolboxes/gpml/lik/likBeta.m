function [varargout] = likBeta(link, hyp, y, mu, s2, inf, i)

% likBeta Beta likelihood function for interval data Y from [0,1].
%
% Report number of hyperparameters
%  s = likBeta ()
%  s = likBeta (link)
%
% Prediction mode
%   lp            = likBeta (link, hyp, y, mu)
%  [lp, ymu, ys2] = likBeta (link, hyp, y, mu, s2)
%
% Inference mode
%  [varargout] = likBeta (link, hyp, y, mu, s2, inf)
%  [varargout] = likBeta (link, hyp, y, mu, s2, inf, i)
%
% Call likFunctions.m to get an explanation of outputs in each mode.
%
% The expression for the likelihood is
%   likBeta(f) = 1/Z * y^(mu*phi-1) * (1-y)^((1-mu)*phi-1) with 
% mean=mu and variance=mu*(1-mu)/(1+phi) where mu = g(f) is the Beta intensity,
% f is a Gaussian process, y is the interval data and
% Z = Gamma(phi)/Gamma(phi*mu)/Gamma(phi*(1-mu)).
% Hence, we have
%   llik(f) = log(likBeta(f)) = -lam*(y-mu)^2/(2*mu^2*y) - log(Zy).
%
% We provide two inverse link functions 'logit' and 'expexp':
%   g(f) = 1/(1+exp(-f)) and g(f) = exp(-exp(-f))).
% The link functions are located at util/glm_invlink_*.m.
%
% Note that for neither link function the likelihood lik(f) is log concave.
% 
% The hyperparameters are:
%
% hyp = [  log(phi)  ]
%
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likFunctions.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme.
%
% Copyright (c) by Hannes Nickisch, 2014-03-04.
%
% See also likFunctions.m.

if nargin<4, varargout = {'1'}; return; end   % report number of hyperparameters

phi = exp(hyp);

if nargin<6                              % prediction mode if inf is not present
  if numel(y)==0,  y = zeros(size(mu)); end
  s2zero = 1; if nargin>4&&numel(s2)>0&&norm(s2)>eps, s2zero = 0; end  % s2==0 ?
  if s2zero                                                    % log probability
    lg = g(mu,link); elg = exp(lg); v = phi*elg; w = phi-v;
    a0 = gammaln(w)-gammaln(phi);
    lp = (v-1).*log(y) + (w-1).*log(1-y) - gammaln(v) - a0;
  else
    lp = likBeta(link, hyp, y, mu, s2, 'infEP');
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
      yv = elg.*(1-elg)/(1+phi);        % second y moment from Beta distribution
      ys2 = (yv+(elg-ymu*oN).^2)*w;
    end
  end
  varargout = {lp,ymu,ys2};
else
  switch inf 
  case 'infLaplace'
    [lg,dlg,d2lg,d3lg] = g(mu,link); elg = exp(lg); v = phi*elg; w = phi-v;
    if nargin<7                                             % no derivative mode
      a0 = gammaln(phi-v)-gammaln(phi);
      lp = (v-1).*log(y) + (w-1).*log(1-y) - gammaln(v) - a0;
      dlp = {}; d2lp = {}; d3lp = {};                         % return arguments
      if nargout>1                           % dlp, derivative of log likelihood
        a1 = v.*(log(y)-log(1-y) + psi(0,w)-psi(0,v));
        dlp = dlg.*a1;
        if nargout>2                    % d2lp, 2nd derivative of log likelihood
          a2 = v.^2.*(psi(1,w)+psi(1,v)); z = dlg.^2+d2lg;
          d2lp = z.*a1 - dlg.^2.*a2;
          if nargout>3                  % d3lp, 3rd derivative of log likelihood
            a3 = v.^3.*(psi(2,w)-psi(2,v));
            d3lp = (dlg.*z+2*dlg.*d2lg+d3lg).*a1 - 3*dlg.*z.*a2 + dlg.^3.*a3;
          end
        end
      end
      varargout = {lp,dlp,d2lp,d3lp};
    else                                                       % derivative mode
                                                  % deriv. of log lik w.r.t. phi
      lp_dhyp = v.*log(y)+w.*log(1-y)-v.*psi(0,v)-w.*psi(0,w)+phi*psi(0,phi);
      a1 = v.*(log(y)-log(1-y) + psi(0,w)-psi(0,v));
      da1 = a1 + v.*(w.*psi(1,w)-v.*psi(1,v));
      dlp_dhyp = dlg.*da1;                                    % first derivative
      a2 = v.^2.*(psi(1,w)+psi(1,v)); z = dlg.^2+d2lg;
      da2 = v.^2.*(w.*psi(2,w)+v.*psi(2,v)) + 2*a2;
      d2lp_dhyp = z.*da1 - dlg.^2.*da2;                      % second derivative
      varargout = {lp_dhyp,dlp_dhyp,d2lp_dhyp};
    end

  case 'infEP'
    if nargin<7                                             % no derivative mode
      % Since we are not aware of an analytical expression of the integral, 
      % we use quadrature.
      varargout = cell(1,nargout);
      [varargout{:}] = lik_epquad({@likBeta,link},hyp,y,mu,s2);
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
  if strcmp(link,'expexp')
    [varargout{:}] = glm_invlink_expexp(f);
  else
    [varargout{:}] = glm_invlink_logit(f);
  end