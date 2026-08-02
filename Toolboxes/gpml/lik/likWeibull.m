function [varargout] = likWeibull(link, hyp, y, mu, s2, inf, i)

% likWeibull - Weibull likelihood function for strictly positive data y. The
% expression for the likelihood is 
%   likWeibull(f) = g1*ka/mu * (g1*y/mu)^(ka-1) * exp(-(g1*y/mu)^ka) with
% gj = gamma(1+j/ka), mean=mu and variance=mu^2*(g2/g1^2-1) where mu = g(f) is
% the Weibull intensity, f is a Gaussian process, y is the positive data.
% Hence, we have llik(f) = log(likWeibull(f)) = 
%                            log(g1*ka/mu) + (ka-1)*log(g1*y/mu) - (g1*y/mu)^ka.
%
% We provide two inverse link functions 'exp' and 'logistic':
%   g(f) = exp(f) and g(f) = log(1+exp(f))).
% The link functions are located at util/glm_invlink_*.m.
%
% Note that for neither link function the likelihood lik(f) is log concave.
% 
% The hyperparameters are:
%
% hyp = [  log(ka)  ]
%
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likFunctions.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme.
%
% Copyright (c) by Hannes Nickisch, 2016-10-04.
%
% See also likFunctions.m.

if nargin<4, varargout = {'1'}; return; end   % report number of hyperparameters

ka = exp(hyp);
lg1 = gammaln(1+1/ka); g1 = exp(lg1); dlg1 = -psi(1+1/ka)/ka;

if nargin<6                              % prediction mode if inf is not present
  if numel(y)==0,  y = zeros(size(mu)); end
  s2zero = 1; if nargin>4&&numel(s2)>0&&norm(s2)>eps, s2zero = 0; end  % s2==0 ?
  if s2zero                                                    % log probability
    lg = g(mu,link);
    lp = lg1 + log(ka) + (ka-1)*(lg1+log(y)) - ka*lg - exp(ka*(lg1+log(y)-lg));
  else
    lp = likWeibull(link, hyp, y, mu, s2, 'infEP');
  end
  ymu = {}; ys2 = {};
  if nargout>1                                 % compute y moments by quadrature
    n = max([length(y),length(mu),length(s2)]); on = ones(n,1);
    N = 20; [t,w] = gauher(N); oN = ones(1,N); lw = ones(n,1)*log(w');
    mu = mu(:).*on; sig = sqrt(s2(:)).*on;                        % vectors only
    lg = g(sig*t'+mu*oN,link); 
    ymu = exp(logsumexp2(lg+lw));     % first moment using Gaussian-Hermite quad
    if nargout>2
      elg = exp(lg); g2 = gamma(1+2/ka);
      yv = elg.^2*(g2/g1^2-1);       % second y moment from Weibull distribution
      ys2 = (yv+(elg-ymu*oN).^2)*w;
    end
  end
  varargout = {lp,ymu,ys2};
else
  switch inf 
  case 'infLaplace'
    [lg,dlg,d2lg,d3lg] = g(mu,link); elg = exp(-ka*lg);
    if nargin<7                                             % no derivative mode
      lp = lg1 + log(ka) + (ka-1)*(lg1+log(y)) -ka*lg - exp(ka*(lg1+log(y)-lg));
      dlp = {}; d2lp = {}; d3lp = {};                         % return arguments
      if nargout>1
        dlp = -ka*dlg + ka*(g1*y).^ka .* elg.*dlg;       % dlp, deriv of log lik
        if nargout>2                    % d2lp, 2nd derivative of log likelihood
          d2lp = -ka*d2lg + ka*(g1*y).^ka .* ( -ka*elg.*dlg.^2 + elg.*d2lg );
          if nargout>3                  % d3lp, 3rd derivative of log likelihood
            a = ka^2*dlg.^3 -3*ka*dlg.*d2lg + d3lg;
            d3lp = - ka*d3lg + ka*(g1*y).^ka .* a .*elg;
          end
        end
      end
      varargout = {lp,dlp,d2lp,d3lp};
    else                                                       % derivative mode
      v = ka*(lg1+log(y)-lg); ev = exp(v);     % derivative of log lik w.r.t. ka
      w = v+ka*dlg1; dw = -ka*dlg; d2w = -ka*d2lg;
      lp_dhyp = 1 + w - ev.*w;
      dlp_dhyp = dw.*(1-ev.*(1+w));                           % first derivative
      d2lp_dhyp = d2w.*(1-ev.*(1+w)) - dw.^2.*(ev.*(2+w));     % and also second
      varargout = {lp_dhyp,dlp_dhyp,d2lp_dhyp};
    end

  case 'infEP'
    if nargin<7                                             % no derivative mode
      % Since we are not aware of an analytical expression of the integral, 
      % we use quadrature.
      varargout = cell(1,nargout);
      [varargout{:}] = lik_epquad({@likWeibull,link},hyp,y,mu,s2);
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