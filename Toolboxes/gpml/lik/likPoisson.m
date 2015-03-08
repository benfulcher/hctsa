function [varargout] = likPoisson(link, hyp, y, mu, s2, inf, i)

% likPoisson - Poisson likelihood function for count data y. The expression for
% the likelihood is 
%   likPoisson(f) = mu^y * exp(-mu) / y! with mean=variance=mu
% where mu = g(f) is the Poisson intensity, f is a
% Gaussian process, y is the non-negative integer count data and 
% y! = gamma(y+1) its factorial. Hence, we have -- with Zy = gamma(y+1) = y! --
%   llik(f) = log(likPoisson(f)) = log(g(f))*y - g(f) - log(Zy).
% The larger the intensity mu, the stronger the likelihood resembles a Gaussian
% since skewness = 1/sqrt(mu) and kurtosis = 1/mu.
%
% We provide two inverse link functions 'exp' and 'logistic':
% For g(f) = exp(f),         we have lik(f) = exp(f*y-exp(f))            / Zy.
% For g(f) = log(1+exp(f))), we have lik(f) = log^y(1+exp(f)))(1+exp(f)) / Zy.
% The link functions are located at util/glm_invlink_*.m.
% 
% Note that for both intensities g(f) the likelihood lik(f) is log concave.
%
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likFunctions.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme.
%
% See also LIKFUNCTIONS.M.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2013-10-16.

if nargin<4, varargout = {'0'}; return; end   % report number of hyperparameters

if nargin<6                              % prediction mode if inf is not present
  if numel(y)==0,  y = zeros(size(mu)); end
  s2zero = 1; if nargin>4, if norm(s2)>0, s2zero = 0; end, end         % s2==0 ?
  if s2zero                                                    % log probability
    lg = g(mu,link);
    lp = lg.*y - exp(lg) - gammaln(y+1);
  else
    lp = likPoisson(link, hyp, y, mu, s2, 'infEP');
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
      yv = elg;                      % second y moment from Poisson distribution
      ys2 = (yv+(elg-ymu*oN).^2)*w;
    end
  end
  varargout = {lp,ymu,ys2};
else
  switch inf 
  case 'infLaplace'
    if nargin<7                                             % no derivative mode
      [lg,dlg,d2lg,d3lg] = g(mu,link); elg = exp(lg);
      lp = lg.*y - elg - gammaln(y+1);
      dlp = {}; d2lp = {}; d3lp = {};                         % return arguments
      if nargout>1
        dlp = dlg.*(y-elg);                  % dlp, derivative of log likelihood
        if nargout>2                    % d2lp, 2nd derivative of log likelihood
          d2lp = d2lg.*(y-elg) - dlg.*dlg.*elg;
          if nargout>3                  % d3lp, 3rd derivative of log likelihood
            d3lp = d3lg.*(y-elg) - dlg.*(dlg.*dlg+3*d2lg).*elg;
          end
        end
      end
      varargout = {lp,dlp,d2lp,d3lp};
    else                                                       % derivative mode
      varargout = {[],[],[]};                         % derivative w.r.t. hypers
    end

  case 'infEP'
    if nargin<7                                             % no derivative mode
      % Since we are not aware of an analytical expression of the integral,
      % hence we use quadrature.
      varargout = cell(1,nargout);
      [varargout{:}] = lik_epquad({@likPoisson,link},hyp,y,mu,s2);
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
  if strcmp(link,'exp')
    [varargout{:}] = glm_invlink_exp(f);
  else
    [varargout{:}] = glm_invlink_logistic(f);
  end