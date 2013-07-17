function [varargout] = likPoisson(kind, hyp, y, mu, s2, inf, i)

% likPoisson - Poisson likelihood function for count data y. The expression for
% the likelihood is 
%   likPoisson(f) = mu^y * exp(-mu) / y! with mean=variance=mu
% where mu = g(f) is the Poisson intensity, f is a
% Gaussian process, y is the non-negative integer count data and 
% y! = gamma(y+1) its factorial. Hence, we have -- with Zy = gamma(y+1) = y! --
%   llik(f) = log(likPoisson(f)) = log(g(f))*y - g(f) - log(Zy).
%
% We provide two kinds of intensities 'exp' and 'logistic':
% For g(f) = exp(f),         we have lik(f) = exp(f*y-exp(f))            / Zy.
% For g(f) = log(1+exp(f))), we have lik(f) = log^y(1+exp(f)))(1+exp(f)) / Zy.
% 
% Note that for both intensities g(f) the likelihood lik(f) is log concave.
%
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likFunctions.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme.
%
% See also LIKFUNCTIONS.M.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2013-01-17

if nargin<4, varargout = {'0'}; return; end   % report number of hyperparameters

if nargin<6                              % prediction mode if inf is not present
  if numel(y)==0,  y = zeros(size(mu)); end
  s2zero = 1; if nargin>4, if norm(s2)>0, s2zero = 0; end, end         % s2==0 ?
  if s2zero                                                    % log probability
    lg = g(mu,kind);
    lp = lg.*y - exp(lg) - gammaln(y+1);
  else
    lp = likPoisson(kind, hyp, y, mu, s2, 'infEP');
  end
  ymu = {}; ys2 = {};
  if nargout>1                                 % compute y moments by quadrature
    n = max([length(y),length(mu),length(s2)]); on = ones(n,1);
    N = 20; [t,w] = gauher(N); oN = ones(1,N); lw = ones(n,1)*log(w');
    mu = mu(:).*on; sig = sqrt(s2(:)).*on;                        % vectors only
    ymu = exp(s(g(sig*t'+mu*oN,kind)+lw));         % exp(mu+s2/2) for kind='exp'
    if nargout>2
      ys2 = ymu;                                  % second y moment equals first
    end
  end
  varargout = {lp,ymu,ys2};
else
  switch inf 
  case 'infLaplace'
    if nargin<7                                             % no derivative mode
      [lg,dlg,d2lg,d3lg] = g(mu,kind); elg = exp(lg);
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
      % we use Gaussian-Hermite quadrature.
      % The following code is generic in the intensity g and relies on the
      % implementation of the infLaplace part.
      n = max([length(y),length(mu),length(s2)]); on = ones(n,1);
      N = 20; [t,w] = gauher(N); oN = ones(1,N); lw = ones(n,1)*log(w');
      y = y(:).*on; mu = mu(:).*on; sig = sqrt(s2(:)).*on;        % vectors only
      [lpi,dlpi,d2lpi] = likPoisson(kind,hyp,y*oN,sig*t'+mu*oN,[],'infLaplace');
      lZ = s(lpi+lw);
      dlZ = {}; d2lZ = {};
      if nargout>1                                     % 1st derivative wrt mean
        % Using p*dlp=dp, p=exp(lp), Z=sum_i wi*pi, dZ = sum_i wi*dpi we obtain
        %   dlZ = sum_i exp(lpi-lZ+lwi)*dlpi = sum_i ai*dlpi.
        a = exp(lpi - lZ*oN + lw);
        dlZ = sum(a.*dlpi,2);
        if nargout>2                                   % 2nd derivative wrt mean
          % Using d2lZ=(d2Z*Z-dZ^2)/Z^2 <=> d2Z=Z*(d2lZ+dlZ^2) and
          % d2Z = sum_i wi*d2Zi, we get d2lZ = sum_i ai*(d2lpi+dlpi^2)-dlZ^2.
          d2lZ = sum(a.*(d2lpi+dlpi.*dlpi),2) - dlZ.*dlZ;
        end
      end
      varargout = {lZ,dlZ,d2lZ};
    else                                                       % derivative mode
      varargout = {[]};                                     % deriv. wrt hyp.lik
    end

  case 'infVB'
    error('infVB not supported')
  end
end

% debug code
% h = 1e-4; ff = [-1e6,-100,-15.01,-14.99,-10,-1,2,100,1e6];
% [lg,dlg,d2lg,d3lg] = g(ff,kind); [lg_h,dlg_h,d2lg_h] = g(ff+h,kind);
% dlgn = (lg_h-lg)/h; d2lgn = (dlg_h-dlg)/h; d3lgn = (d2lg_h-d2lg)/h;
% 
% [dlg(:),dlgn(:), d2lg(:),d2lgn(:), d3lg(:),d3lgn(:)]
% keyboard

% compute the log Poisson intensity
function varargout = g(f,kind)
  varargout = cell(nargout, 1);  % allocate the right number of output arguments
  if strcmp(kind,'exp')
    [varargout{:}] = g_exp(f);
  else
    [varargout{:}] = g_logistic(f);
  end

% compute the log Poisson intensity for g(f) = exp(f)
function [lg,dlg,d2lg,d3lg] = g_exp(f)
  lg = f;
  if nargout>1
    dlg = ones(size(f));
    if nargout>2
      d2lg = zeros(size(f));
      if nargout>2
        d3lg = zeros(size(f));
      end
    end
  end

% compute the log Poisson intensity for g(f) = log(1+exp(f)))
function [lg,dlg,d2lg,d3lg] = g_logistic(f)
  l1pef = max(0,f) + log(1+exp(-abs(f)));         % safely compute log(1+exp(f))
  lg = log(l1pef); id = f<-15; lg(id) = f(id);   % fix log(log(1+exp(f))) limits
  if nargout>1
    sm = 1./(1+exp(-f));
    dlg = sm./l1pef; dlg(f<-15) = 1;
    if nargout>2
      sp = 1./(1+exp(f));
      d2lg = dlg.*(sp-dlg);
      if nargout>2
        d3lg = d2lg.*(sp-2*dlg) - dlg.*sp.*sm;
      end
    end
  end

% computes y = log( sum(exp(x),2) ), the softmax in a numerically safe way by
%  subtracting the row maximum to avoid cancelation after taking the exp
%  the sum is done along the rows
function [y,x] = s(logx)
  N = size(logx,2); max_logx = max(logx,[],2);
  % we have all values in the log domain, and want to calculate a sum
  x = exp(logx-max_logx*ones(1,N));
  y = log(sum(x,2)) + max_logx;
