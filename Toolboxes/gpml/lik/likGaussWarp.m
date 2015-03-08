function [varargout] = likGaussWarp(warp, hyp, y, mu, varargin)

% likGaussWarp - Warped Gaussian likelihood for regression. 
% The expression for the likelihood is 
%   likGaussWarp( y | t ) = likGauss( g(y) | t ) * g'(y),
% where likGauss is the Gaussian likelihood and g is the warping function.
%
% The hyperparameters are:
%
% hyp = [ theta_1
%         theta_2
%           ..
%         theta_ng
%         log(sn) ]
%
% Here, sn is the standard deviation of the underlying Gaussian and theta_i for
% i=1..ng are the ng hyperparameters of the warping function g.
%
% At the moment, likGaussWarp offers 3 different warping functions:
% id                   yields g(y) = y  =>  likGaussWarp = likGauss
% poly<m> e.g. 'poly1' yields g(y) = y  =>  likGaussWarp = likGauss
%              'poly3' yields g(y) = y + c1*sy*ay^2 + c2*sy*ay^3
%                             where sy = sign(y), ay = abs(y), cj = exp(theta_j)
% tanh<m> e.g. 'tanh0' yields g(y) = y  =>  likGaussWarp = likGauss
%              'tanh2' yields g(y) = y + a1*tanh(b1*(y+c1)) + a2*tanh(b2*(y+c2))
%                 where aj = exp(theta_j), bj = exp(theta_j+m), bj = theta_j+2*m
%
% The code is based on the exposition in the paper Warped Gaussian Processes,
% NIPS, 2003 by Edward Snelson, Carl Rasmussen and Zoubin Ghahramani.
%
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likFunctions.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme.
%
% Copyright (c) by Hannes Nickisch, 2013-10-24.
%
% See also LIKFUNCTIONS.M.

lik = {@likGauss}; % in principle any likelihood function can be warped but only
% for homoscedastic likelihoods, in particular Gaussian has feasible integrals
if numel(warp)==0, warp = 'id'; end               % set default warping function
ng = g(warp);               % number of hyperparameters for the warping function
nhyp = ['(',num2str(ng),'+',feval(lik{:}),')'];      % number of hyperparameters
if nargin<4, varargout = {nhyp}; return, end       % report number of parameters
nhyp = eval(nhyp);
if nhyp>length(hyp), error('not enough hyperparameters'), end

[gy,lgpy] = g(warp,y,hyp(1:ng));                     % evaluate warping function
i = 0; if nargin>6, i = varargin{3}; varargin{3} = varargin{3}-ng; end
varargout = cell(nargout,1);              % allocate memory for output arguments
if i==0 || ng<i                               % only evaluate the required parts
  [varargout{:}] = feval(lik{:},hyp(ng+1:end),gy,mu,varargin{:});     % eval lik
end
if nargin<6                              % prediction mode if inf is not present
  if numel(y)==0,  y = zeros(size(mu)); end
  s2zero = 1;                                                          % s2==0 ?
  if nargin>4, s2 = varargin{1}; if norm(s2)>0, s2zero = 0; end, end
  if s2zero                                                    % log probability
    lp = likGaussWarp(warp, hyp, y, mu, [], 'infLaplace'); s2 = 0*mu;
  else
    lp = likGaussWarp(warp, hyp, y, mu, s2, 'infEP');               % prediction
  end
  if nargout>0, varargout{1} = lp; end                       % works for any lik
  % the predictive moments are very hard to compute for lik not being likGauss
  if nargout>1
    ymu = mu;                                                % first g(y) moment
    sn2 = exp(2*hyp(ng+1));                            % Gaussian noise variance
    ys2 = s2 + sn2;                                         % second g(y) moment
%     ymuM = ig(warp,ymu,hyp(1:ng));                                    % median
%     yupp = ig(warp,ymu+2*sqrt(ys2),hyp(1:ng));       % 95% confidence interval
%     ylow = ig(warp,ymu-2*sqrt(ys2),hyp(1:ng));
%     ys2C = (yupp-ylow).^2/16;
    N = 20; [t,w] = gauher(N); oN = ones(1,N);     % Gaussian-Hermite quadrature
    Z = sqrt(ys2(:))*t'+ymu(:)*oN;
    Y = ig(warp,Z,hyp(1:ng));
    ymu = Y*w; ys2 = (Y-ymu*oN).^2*w;                % first and second y moment
    varargout{2} = reshape(ymu,size(mu));
    if nargout>2
      varargout{3} = ys2;
    end
  end

else
  inf = varargin{2};                                    % obtain input variables
  switch inf
  case {'infLaplace','infEP'}                     % they have the same structure
    if nargin<7                                             % no derivative mode
      if nargout>0, varargout{1} = varargout{1} + lgpy; end
    else                                                       % derivative mode
      if i<=ng                  % derivatives w.r.t. warping function parameters
        n = max([numel(y),numel(mu)]);
        for j=2:nargout, varargout{j} = zeros(n,1); end
        [dgy,dlgpy] = g(warp,y,hyp(1:ng),i);       % warping function derivative
        out = cell(nargout+1,1);                               % allocate memory
        [out{:}] = likGaussWarp(warp, hyp, y, mu, varargin{1:2});    % query lik
        % works only for homoscedastic likelihoods where y and mu can be swapped
        if nargout>0, varargout{1} = dlgpy - out{2}.*dgy; end % apply chain rule
        if nargout>1, varargout{2} =       - out{3}.*dgy; end
        if nargout>2, varargout{3} =       - out{4}.*dgy; end
      end
    end

  case 'infVB'           % output does not depend on mu and following parameters
  end
end

% marshalling of parameters and available warping functions
function varargout = g(warp,varargin)
  varargout = cell(nargout, 1);  % allocate the right number of output arguments
  if strcmp(warp,'id')               % indentity warping likGaussWarp = likGauss
    if nargin<2
      if nargout>0, varargout{1} = 0; end
    elseif nargin<4
      if nargout>0, varargout{1} =   varargin{1}; end
      if nargout>1, varargout{2} = 0*varargin{1}; end
    end
  elseif numel(strfind(warp,'poly'))>0
    m = str2double(warp(5:end));
    if nargin<2 && nargout>0, varargout{1} = m-1; return, end
    [varargout{:}] = g_poly(varargin{:});
  elseif numel(strfind(warp,'tanh'))>0
    m = str2double(warp(5:end));
    if nargin<2 && nargout>0, varargout{1} = 3*m; return, end
    [varargout{:}] = g_tanh(varargin{:});
  end

% invert g(y) = z <=> ig(z) = y via bisection search + Newton iterations
function [y,n,d] = ig(warp,z,hyp)
  y = z; gy = g(warp,z,hyp)-z; dz = max(abs(z(:)));    % lower bound search ylow
  while any(0<gy(:)), y(0<gy) = y(0<gy)-dz; gy = g(warp,y,hyp)-z; end, ylow = y;
  y = z; gy = g(warp,z,hyp)-z; dz = max(abs(z(:)));    % upper bound search yupp
  while any(0>gy(:)), y(0>gy) = y(0>gy)+dz; gy = g(warp,y,hyp)-z; end, yupp = y;
  for n=1:12                                    % bisection search ylow<=y<=yupp
    d = max(abs(gy(:))); if d<sqrt(eps), break, end
    y = (ylow+yupp)/2; gy = g(warp,y,hyp)-z;
    ylow(gy<0) = y(gy<0); yupp(gy>0) = y(gy>0);
  end
  for n=1:12                                                 % Newton iterations
    [gy,lgpy] = g(warp,y,hyp); gpy = exp(lgpy);
    y = y - (gy-z)./gpy;
    y(y<ylow) = ylow(y<ylow); y(y>yupp) = yupp(y>yupp);          % keep brackets
    d = max( abs(gy(:)-z(:)) );
    if d<sqrt(eps), break, end
  end
  if n==10 || d>sqrt(eps), fprintf('Not converged: res=%1.4e\n',d), end

% poly warping function g(y) and log of the derivative log(g'(y))>0
% or derivatives of the latter w.r.t. ith hyperparameter
function [gy,lgpy] = g_poly(y,hyp,i)
  m = numel(hyp)+1;
  c = exp(hyp);
  if nargin==2                                                 % function values
    gy = y; gpy = 1; ay = abs(y);
    for j=2:m
      gy  = gy + c(j-1)*ay.^j;
      gpy = gpy + c(j-1)*j*ay.^(j-1);
    end
    gy = sign(y).*gy;
    lgpy = log(gpy);
  else                                                             % derivatives
    gpy = 1; ay = abs(y);
    for j=2:m
      gpy = gpy + c(j-1)*j*ay.^(j-1);
    end
    gy = c(i)*ay.^j;
    lgpy = c(i)*j*ay.^(j-1)./gpy;
  end

% tanh warping function g(y) and log of the derivative log(g'(y))>0
% or derivatives of the latter w.r.t. ith hyperparameter
function [gy,lgpy] = g_tanh(y,hyp,i)
  m = numel(hyp)/3;
  a = exp(hyp(1:m)); b = exp(hyp(m+(1:m))); c = hyp(2*m+(1:m));
  if nargin==2                                                 % function values
    gy = y; gpy = 1;
    for j=1:m
      ai = a(j); bi = b(j); ci = c(j); ti = tanh(bi*(y+ci)); dti = 1-ti.^2;
      gy  = gy  + ai    *ti;
      gpy = gpy + ai*bi*dti;
    end
    lgpy = log(gpy);
  else                                                             % derivatives
    gpy = 1;
    for j=1:m
      ai = a(j); bi = b(j); ci = c(j); ti = tanh(bi*(y+ci)); dti = 1-ti.^2;
      gpy = gpy + ai*bi*dti;
    end
    if i<=m
      j = i;
      ai = a(j); bi = b(j); ci = c(j); ti = tanh(bi*(y+ci)); dti = 1-ti.^2;
      gy = ai*ti;
      lgpy = ai*bi*dti./gpy;
    elseif i<=2*m
      j = i-m;
      ai = a(j); bi = b(j); ci = c(j); ti = tanh(bi*(y+ci)); dti = 1-ti.^2;
      gy = ai*bi*dti.*(y+ci);
      d2ti = -2*ti.*dti;
      lgpy = ai*bi*(dti+bi*d2ti.*(y+ci))./gpy;
    else
      j = i-2*m;
      ai = a(j); bi = b(j); ci = c(j); ti = tanh(bi*(y+ci)); dti = 1-ti.^2;
      gy = ai*bi*dti;
      d2ti = -2*ti.*dti;
      lgpy = ai*bi^2*d2ti./gpy;
    end
  end