function [varargout] = likLogGauss(link, hyp, y, mu, s2, inf, i)

% likLogGauss - Log Gaussian likelihood function for strictly positive data y.
% The expression for the likelihood is 
%   likLogGauss(f) = 1/Zy*exp( -[log(mu)-sn^2/2 - log(y)]^2 / (2*sn^2) ) with 
% mean=mu and variance=mu^3/lam where mu = g(f) is the Log Gaussian 
% intensity, f is a Gaussian process, y is the strictly positive data and
% Zy = y*sn*sqrt(2*pi) is the normalizer. Note that, we have
% log(y) ~ N(mu-sn^2/2,sn^2)
%
% Hence, we have
%   llik(f) = log(likLogGauss(f)) = -(log(mu)-sn^2/2-log(y))^2/(2*sn^2)-log(Zy).
%
% We provide two inverse link functions 'exp' and 'logistic':
%   g(f) = exp(f) and g(f) = log(1+exp(f))).
% The link functions are located at util/glm_invlink_*.m.
%
% Note that the likelihood lik(f) is log concave for the exponential link.
%
% The hyperparameters are:
%
% hyp = [  log(sn)  ]
%
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likFunctions.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme.
%
% Copyright (c) by Hannes Nickisch, 2018-09-19.
%
% See also likFunctions.m.

if nargin<4, varargout = {'1'}; return; end   % report number of hyperparameters

sn = exp(hyp);

if nargin<6                              % prediction mode if inf is not present
  if numel(y)==0,  y = zeros(size(mu)); end
  s2zero = 1; if nargin>4&&numel(s2)>0&&norm(s2)>eps, s2zero = 0; end  % s2==0 ?
  if s2zero                                                    % log probability
    lg = g(mu,link);
    lZy = log(y*sn*sqrt(2*pi));                         % normalisation constant
    lp = -(lg-sn^2/2 - log(y)).^2 / (2*sn^2) - lZy;
  else
    lp = likLogGauss(link, hyp, y, mu, s2, 'infEP');
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
      yv = (exp(sn^2)-1)*elg.^2;% second y moment from log Gaussian distribution
      ys2 = (yv+(elg-ymu*oN).^2)*w;
    end
  end
  varargout = {lp,ymu,ys2};
else
  switch inf 
  case 'infLaplace'
    [lg,dlg,d2lg,d3lg] = g(mu,link); elg = exp(lg);
    if nargin<7                                             % no derivative mode
      lZy = log(y*sn*sqrt(2*pi));                       % normalisation constant
      a = lg-sn^2/2 - log(y); lp = -a.^2 / (2*sn^2) - lZy;
      dlp = {}; d2lp = {}; d3lp = {};                         % return arguments
      if nargout>1
        dlp = -2*dlg.*a / (2*sn^2);   % dlp, deriv of log lik
        if nargout>2                    % d2lp, 2nd derivative of log likelihood
          d2lp = -2*(d2lg.*a+dlg.^2) / (2*sn^2);
          if nargout>3                  % d3lp, 3rd derivative of log likelihood
            d3lp = -2*(d3lg.*a+2*dlg.*d2lg) / (2*sn^2);
          end
        end
      end
      varargout = {lp,dlp,d2lp,d3lp};
    else                                                       % derivative mode
      a = lg-sn^2/2 - log(y); da = -sn^2;
      lp_dhyp = a.*(a-da) / sn^2 - 1;          % derivative of log lik w.r.t. sn
      dlp_dhyp = dlg.*(2*a-da) / sn^2;                        % first derivative
      d2lp_dhyp = d2lg.*(2*a-da) / sn^2 + 2*dlg.*dlg / sn^2;    % 2nd derivative
      varargout = {lp_dhyp,dlp_dhyp,d2lp_dhyp};
    end

  case 'infEP'
    if nargin<7                                             % no derivative mode
      % Since we are not aware of an analytical expression of the integral, 
      % we use quadrature.
      varargout = cell(1,nargout);
      [varargout{:}] = lik_epquad({@likLogGauss,link},hyp,y,mu,s2);
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