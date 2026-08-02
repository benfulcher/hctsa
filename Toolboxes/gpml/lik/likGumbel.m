function [varargout] = likGumbel(sign, hyp, y, mu, s2, inf, i)

% likGumbel - Gumbel likelihood function for extremal value regression. 
% The expression for the likelihood is
%   likGumbel(t) = exp(-z-exp(-z))/be, z = ga+s*(y-t)/be, be = sn*sqrt(6)/pi
% where s={+1,-1} is a sign switching between left and right skewed, ga is the
% Euler-Mascheroni constant, y is the mean, sn^2 is the variance.
% The skewness and kurtosis of likGumbel are 1.14*s and 2.4, respectively. 
%
% The hyperparameters are:
%
% hyp = [ log(sn)  ]
%
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likFunctions.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme.
%
% Copyright (c) by Hannes Nickisch, 2013-11-01.
%
% See also likFunctions.m.

if nargin<4, varargout = {'1'}; return; end   % report number of hyperparameters
if sign=='-', s = -1; else s = 1; end                 % extract sign of skewness

sn2 = exp(2*hyp);                                      % extract hyperparameters
ga = 0.5772156649;                                   % Euler-Mascheroni constant
be = sqrt(6*sn2)/pi;
lZ = -log(be);

if nargin<6                              % prediction mode if inf is not present
  if numel(y)==0,  y = zeros(size(mu)); end
  s2zero = 1; if nargin>4&&numel(s2)>0&&norm(s2)>eps, s2zero = 0; end  % s2==0 ?
  if s2zero                                         % log probability evaluation
    lp = likGumbel(sign, hyp, y, mu, [], 'infLaplace'); s2 = 0;
  else                                                              % prediction
    lp = likGumbel(sign, hyp, y, mu, s2, 'infEP');
  end
  ymu = {}; ys2 = {};
  if nargout>1
    ymu = mu;                                                   % first y moment
    if nargout>2
      ys2 = s2 + sn2;                                          % second y moment
    end
  end
  varargout = {lp,ymu,ys2};
else
  switch inf 
  case 'infLaplace'
    z = ga+s*(y-mu)/be; emz = exp(-z);
    if nargin<7                                             % no derivative mode
      dlp = {}; d2lp = {}; d3lp = {};
      lp = lZ -z -emz;
      if nargout>1
        dz = -s/be;                                                     % dz/dmu
        dlp = dz*(emz-1);                    % dlp, derivative of log likelihood
        if nargout>2                    % d2lp, 2nd derivative of log likelihood
          d2lp = -dz^2*emz;
          if nargout>3                  % d3lp, 3rd derivative of log likelihood
            d3lp = dz^3*emz;
          end
        end
      end
      varargout = {lp,dlp,d2lp,d3lp};
    else                                             % derivative w.r.t. log(sn)
      dz = -s/be;                                                       % dz/dmu
      dzs = -s*(y-mu)/be;                                          % dz/dlog(sn)
      lp_dhyp   =  dzs.*(emz-1) -1;
      dlp_dhyp  = dz*(1-emz.*(1+dzs));
      d2lp_dhyp = dz^2*emz.*(2+dzs);
      varargout = {lp_dhyp,dlp_dhyp,d2lp_dhyp};
    end

  case 'infEP'
    if nargout>1
      error('infEP not supported since likT is not log-concave')
    end
    n = max([length(y),length(mu),length(s2)]); on = ones(n,1);
    y = y(:).*on; mu = mu(:).*on; sig = sqrt(s2(:)).*on;          % vectors only
    % since we are not aware of an analytical expression of the integral, 
    % we use Gaussian-Hermite quadrature
    N = 20; [t,w] = gauher(N); oN = ones(1,N);
    lZ = likGumbel(sign, hyp, y*oN, sig*t'+mu*oN, []);
    lZ = log_expA_x(lZ,w); % log( exp(lZ)*w )
    varargout = {lZ};

  case 'infVB'
    error('infVB not supported')
  end
end

%  computes y = log( exp(A)*x ) in a numerically safe way by subtracting the
%  maximal value in each row to avoid cancelation after taking the exp
function y = log_expA_x(A,x)
  N = size(A,2);  maxA = max(A,[],2);      % number of columns, max over columns
  y = log(exp(A-maxA*ones(1,N))*x) + maxA;  % exp(A) = exp(A-max(A))*exp(max(A))