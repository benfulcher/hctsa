function [varargout] = likErf(hyp, y, mu, s2, inf, i)

% likErf Error function or cumulative Gaussian likelihood function for binary
% classification or probit regression
%
% Report number of hyperparameters
%  s = likErf ()
%  s = likErf (link)
%
% Prediction mode
%   lp            = likErf (hyp, y, mu)
%  [lp, ymu, ys2] = likErf (hyp, y, mu, s2)
%
% Inference mode
%  [varargout] = likErf (hyp, y, mu, s2, inf)
%  [varargout] = likErf (hyp, y, mu, s2, inf, i)
%
% Call likFunctions.m to get an explanation of outputs in each mode.
%
% The expression for the likelihood is 
%   likErf(t) = (1+erf(t/sqrt(2)))/2 = normcdf(t).
%
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likFunctions.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme.
% 
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2014-03-19.
%
% See also likFunctions.m.

if nargin<3, varargout = {'0'}; return; end   % report number of hyperparameters
if nargin>1, y = sign(y); y(y==0) = 1; else y = 1; end % allow only +/- 1 values
if numel(y)==0, y = 1; end

if nargin<5                              % prediction mode if inf is not present
  y = y.*ones(size(mu));                                       % make y a vector
  s2zero = 1; if nargin>3&&numel(s2)>0&&norm(s2)>eps, s2zero = 0; end  % s2==0 ?
  if s2zero                                         % log probability evaluation
    lp = logphi(y.*mu);
  else                                                              % prediction
    lp = likErf(hyp, y, mu, s2, 'infEP');
  end
  p = exp(lp); ymu = {}; ys2 = {};
  if nargout>1
    ymu = 2*p-1;                                                % first y moment
    if nargout>2
      ys2 = 4*p.*(1-p);                                        % second y moment
    end
  end
  varargout = {lp,ymu,ys2};
else                                                            % inference mode
  switch inf 
  case 'infLaplace'
    if nargin<6                                             % no derivative mode
      f = mu; yf = y.*f;                            % product latents and labels
      varargout = cell(nargout,1); [varargout{:}] = logphi(yf);   % query logphi
      if nargout>1
        varargout{2} = y.*varargout{2};
        if nargout>3, varargout{4} = y.*varargout{4}; end
      end
    else                                                       % derivative mode
      varargout = {[],[],[]};                         % derivative w.r.t. hypers
    end

  case 'infEP'
    if nargin<6                                             % no derivative mode
      z = mu./sqrt(1+s2); dlZ = {}; d2lZ = {};
      if numel(y)>0, z = z.*y; end
      if nargout<=1, lZ = logphi(z);                         % log part function
      else          [lZ,n_p] = logphi(z); end
      if nargout>1
        if numel(y)==0, y=1; end
        dlZ = y.*n_p./sqrt(1+s2);                      % 1st derivative wrt mean
        if nargout>2, d2lZ = -n_p.*(z+n_p)./(1+s2); end         % 2nd derivative
      end
      varargout = {lZ,dlZ,d2lZ};
    else                                                       % derivative mode
      varargout = {[]};                                     % deriv. wrt hyp.lik
    end

  case 'infVB'
    error('infVB not supported')
  end
end