function [varargout] = likT(hyp, y, mu, s2, inf, i)

% likT - Student's t likelihood function for regression. 
% The expression for the likelihood is
%   likT(t) = Z * ( 1 + (t-y)^2/(nu*sn^2) ).^(-(nu+1)/2),
% where Z = gamma((nu+1)/2) / (gamma(nu/2)*sqrt(nu*pi)*sn)
% and y is the mean (for nu>1) and nu*sn^2/(nu-2) is the variance (for nu>2).
%
% The hyperparameters are:
%
% hyp = [ log(nu-1)
%         log(sn)  ]
%
% Note that the parametrisation guarantees nu>1, thus the mean always exists.
%
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likFunctions.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2013-09-02.
%
% See also LIKFUNCTIONS.M.

if nargin<3, varargout = {'2'}; return; end   % report number of hyperparameters

numin = 1;                                                 % minimum value of nu
nu = exp(hyp(1))+numin; sn2 = exp(2*hyp(2));           % extract hyperparameters
lZ = gammaln(nu/2+1/2) - gammaln(nu/2) - log(nu*pi*sn2)/2;

if nargin<5                              % prediction mode if inf is not present
  if numel(y)==0,  y = zeros(size(mu)); end
  s2zero = 1; if nargin>3, if norm(s2)>0, s2zero = 0; end, end         % s2==0 ?
  if s2zero                                         % log probability evaluation
    lp = lZ - (nu+1)*log( 1+(y-mu).^2./(nu.*sn2) )/2; s2 = 0;
  else                                                              % prediction
    lp = likT(hyp, y, mu, s2, 'infEP');
  end
  ymu = {}; ys2 = {};
  if nargout>1
    ymu = mu;                       % first y moment; for nu<=1 this is the mode
    if nargout>2
      if nu<=2
        ys2 = Inf(size(mu));                    % variance does not always exist
      else
        ys2 = (s2 + nu*sn2/(nu-2)).*ones(size(mu));            % second y moment
      end
    end
  end
  varargout = {lp,ymu,ys2};
else
  switch inf 
  case 'infLaplace'
    r = y-mu; r2 = r.*r;
    if nargin<6                                             % no derivative mode
      dlp = {}; d2lp = {}; d3lp = {};
      lp = lZ - (nu+1)*log( 1+r2./(nu.*sn2) )/2;
      if nargout>1
        a = r2+nu*sn2;
        dlp = (nu+1)*r./a;                   % dlp, derivative of log likelihood
        if nargout>2                    % d2lp, 2nd derivative of log likelihood
          d2lp = (nu+1)*(r2-nu*sn2)./a.^2;
          if nargout>3                  % d3lp, 3rd derivative of log likelihood
            d3lp = (nu+1)*2*r.*(r2-3*nu*sn2)./a.^3;
          end
        end
      end
      varargout = {lp,dlp,d2lp,d3lp};
    else                                                       % derivative mode
      a = r2+nu*sn2; a2 = a.*a; a3 = a2.*a;
      if i==1                                             % derivative w.r.t. nu
        lp_dhyp =  nu*( psi(nu/2+1/2)-psi(nu/2) )/2 - 1/2 ...
                  -nu*log(1+r2/(nu*sn2))/2 +(nu/2+1/2)*r2./(nu*sn2+r2);
        lp_dhyp = (1-numin/nu)*lp_dhyp;          % correct for lower bound on nu
        dlp_dhyp = nu*r.*( a - sn2*(nu+1) )./a2;
        dlp_dhyp = (1-numin/nu)*dlp_dhyp;        % correct for lower bound on nu       
        d2lp_dhyp = nu*( r2.*(r2-3*sn2*(1+nu)) + nu*sn2^2 )./a3;
        d2lp_dhyp = (1-numin/nu)*d2lp_dhyp;      % correct for lower bound on nu
      else                                                % derivative w.r.t. sn
        lp_dhyp   =  (nu+1)*r2./a - 1;
        dlp_dhyp  = -(nu+1)*2*nu*sn2*r./a2;
        d2lp_dhyp =  (nu+1)*2*nu*sn2*(a-4*r2)./a3;
      end
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
    lZ = likT(hyp, y*oN, sig*t'+mu*oN, []);
    lZ = log_expA_x(lZ,w); % log( exp(lZ)*w )
    varargout = {lZ};

  case 'infVB'
    % variational lower site bound
    % t(s) \propto (1+(s-y)^2/(nu*s2))^(-nu/2+1/2)
    % the bound has the form: (b+z/ga)*f - f.^2/(2*ga) - h(ga)/2
    n = numel(s2); b = zeros(n,1); y = y.*ones(n,1); z = y;
    varargout = {b,z};
  end
end

%  computes y = log( exp(A)*x ) in a numerically safe way by subtracting the
%  maximal value in each row to avoid cancelation after taking the exp
function y = log_expA_x(A,x)
  N = size(A,2);  maxA = max(A,[],2);      % number of columns, max over columns
  y = log(exp(A-maxA*ones(1,N))*x) + maxA;  % exp(A) = exp(A-max(A))*exp(max(A))
