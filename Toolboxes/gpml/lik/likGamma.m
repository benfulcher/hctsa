function [varargout] = likGamma(link, hyp, y, mu, s2, inf, i)

% likGamma - Gamma likelihood function for strictly positive data y. The
% expression for the likelihood is 
%   likGamma(f) = al^al*y^(al-1)/gamma(al) * exp(-y*al/mu) / mu^al with 
% mean=mu and variance=mu^2/al where mu = g(f) is the Gamma intensity, f is a
% Gaussian process, y is the strictly positive data. Hence, we have -- with
% log(Zy) = log(gamma(al)) - al*log(al) + (1-al)*log(y)
%   llik(f) = log(likGamma(f)) = -al*( log(g(f)) + y/g(f) ) - log(Zy).
% The larger one chooses al, the stronger the likelihood resembles a Gaussian
% since skewness = 2/sqrt(al) and kurtosis = 6/al.
%
% We provide two inverse link functions 'exp' and 'logistic':
%   g(f) = exp(f) and g(f) = log(1+exp(f))).
% The link functions are located at util/glm_invlink_*.m.
%
% Note that for neither link function the likelihood lik(f) is log concave.
% 
% The hyperparameters are:
%
% hyp = [  log(al)  ]
%
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likFunctions.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme.
%
% Copyright (c) by Hannes Nickisch, 2016-10-04.
%
% See also likFunctions.m.

if nargin<4, varargout = {'1'}; return; end   % report number of hyperparameters

al = exp(hyp);

if nargin<6                              % prediction mode if inf is not present
  if numel(y)==0,  y = zeros(size(mu)); end
  s2zero = 1; if nargin>4&&numel(s2)>0&&norm(s2)>eps, s2zero = 0; end  % s2==0 ?
  if s2zero                                                    % log probability
    lg = g(mu,link);
    lZy = gammaln(al) - al*log(al) + (1-al)*log(y);     % normalisation constant
    lp = -al*(lg+y./exp(lg)) - lZy;
  else
    lp = likGamma(link, hyp, y, mu, s2, 'infEP');
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
      yv = elg.^2/al;                  % second y moment from Gamma distribution
      ys2 = (yv+(elg-ymu*oN).^2)*w;
    end
  end
  varargout = {lp,ymu,ys2};
else
  switch inf 
  case 'infLaplace'
    [lg,dlg,d2lg,d3lg] = g(mu,link); elg = exp(lg);
    if nargin<7                                             % no derivative mode
      lZy = gammaln(al) - al*log(al) + (1-al)*log(y);   % normalisation constant
      lp = -al*(lg+y./elg) - lZy;
      dlp = {}; d2lp = {}; d3lp = {};                         % return arguments
      if nargout>1
        dlp = -al*dlg.*(1-y./elg);           % dlp, derivative of log likelihood
        if nargout>2                    % d2lp, 2nd derivative of log likelihood
          d2lp = -al*d2lg.*(1-y./elg) - al*dlg.*dlg.*y./elg;
          if nargout>3                  % d3lp, 3rd derivative of log likelihood
            d3lp = -al*d3lg.*(1-y./elg) + al*dlg.*(dlg.*dlg-3*d2lg).*y./elg;
          end
        end
      end
      varargout = {lp,dlp,d2lp,d3lp};
    else                                                       % derivative mode
      dlZy = al*psi(0,al) - al*(log(al) + 1 + log(y));
      lp_dhyp = -al*(lg+y./elg) - dlZy; % derivative of log likelihood w.r.t. al
      dlp_dhyp = -al*dlg.*(1-y./elg);                         % first derivative
      d2lp_dhyp = -al*d2lg.*(1-y./elg) - al*dlg.*dlg.*y./elg;  % and also second
      varargout = {lp_dhyp,dlp_dhyp,d2lp_dhyp};
    end

  case 'infEP'
    if nargin<7                                             % no derivative mode
      % Since we are not aware of an analytical expression of the integral, 
      % we use quadrature.
      varargout = cell(1,nargout);
      [varargout{:}] = lik_epquad({@likGamma,link},hyp,y,mu,s2);
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