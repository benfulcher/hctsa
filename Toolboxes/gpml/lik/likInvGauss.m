function [varargout] = likInvGauss(link, hyp, y, mu, s2, inf, i)

% likInvGauss - Inverse Gaussian likelihood function for strictly positive data
% y. The expression for the likelihood is 
%   likInvGauss(f) = sqrt(lam/(2*pi*y^3))*exp(-lam*(mu-y)^2/(2*mu^2*y)) with 
% mean=mu and variance=mu^3/lam where mu = g(f) is the Inverse Gaussian 
% intensity, f is a Gaussian process, y is the strictly positive data. 
% Hence, we have -- with log(Zy) = -(log(lam)-log(2*pi*y^3))/2
%   llik(f) = log(likInvGauss(f)) = -lam*(y-mu)^2/(2*mu^2*y) - log(Zy).
% The larger one chooses lam, the stronger the likelihood resembles a Gaussian
% since skewness = 3*sqrt(mu/lam) and kurtosis = 15*mu/lam.
%
% We provide two inverse link functions 'exp' and 'logistic':
%   g(f) = exp(f) and g(f) = log(1+exp(f))).
% The link functions are located at util/glm_invlink_*.m.
%
% Note that for neither link function the likelihood lik(f) is log concave.
%
% The hyperparameters are:
%
% hyp = [  log(lam)  ]
%
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likFunctions.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme.
%
% See also LIKFUNCTIONS.M.
%
% Copyright (c) by Hannes Nickisch, 2013-10-16.

if nargin<4, varargout = {'1'}; return; end   % report number of hyperparameters

lam = exp(hyp);

if nargin<6                              % prediction mode if inf is not present
  if numel(y)==0,  y = zeros(size(mu)); end
  s2zero = 1; if nargin>4, if norm(s2)>0, s2zero = 0; end, end         % s2==0 ?
  if s2zero                                                    % log probability
    lg = g(mu,link); elg = exp(lg);
    lZy = -(log(lam)-log(2*pi*y.^3))/2;                 % normalisation constant
    lp = -lam*(y./elg-1).^2 ./(2*y) - lZy;
  else
    lp = likInvGauss(link, hyp, y, mu, s2, 'infEP');
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
      yv = elg.^3/lam;      % second y moment from inverse Gaussian distribution
      ys2 = (yv+(elg-ymu*oN).^2)*w;
    end
  end
  varargout = {lp,ymu,ys2};
else
  switch inf 
  case 'infLaplace'
    [lg,dlg,d2lg,d3lg] = g(mu,link); elg = exp(lg);
    if nargin<7                                             % no derivative mode
      lZy = -(log(lam)-log(2*pi*y.^3))/2;               % normalisation constant
      lp = -lam*(y./elg-1).^2 ./(2*y) - lZy;
      dlp = {}; d2lp = {}; d3lp = {};                         % return arguments
      if nargout>1
        dlp = lam*(y./elg-1).*dlg./elg;      % dlp, derivative of log likelihood
        if nargout>2                    % d2lp, 2nd derivative of log likelihood
          d2lp = lam*( (y./elg-1).*(d2lg-dlg.^2) - y.*dlg.^2./elg )./elg;
          if nargout>3                  % d3lp, 3rd derivative of log likelihood
            d3lp = lam*(  (y./elg-1) .* (4*dlg.^3-6*dlg.*d2lg+d3lg)  ...
                                      + (3*dlg.^3 - 3*dlg.*d2lg)        )./ elg;
          end
        end
      end
      varargout = {lp,dlp,d2lp,d3lp};
    else                                                       % derivative mode
      lp_dhyp = 1/2-lam*(y./elg-1).^2 ./(2*y);    % deriv. of log lik w.r.t. lam
      dlp_dhyp = lam*(y./elg-1).*dlg./elg;                    % first derivative
      d2lp_dhyp = lam*( (y./elg-1).*(d2lg-dlg.^2) - y.*dlg.^2./elg )./elg; % 2nd
      varargout = {lp_dhyp,dlp_dhyp,d2lp_dhyp};
    end

  case 'infEP'
    if nargin<7                                             % no derivative mode
      % Since we are not aware of an analytical expression of the integral, 
      % we use quadrature.
      varargout = cell(1,nargout);
      [varargout{:}] = lik_epquad({@likInvGauss,link},hyp,y,mu,s2);
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