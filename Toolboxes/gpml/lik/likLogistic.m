function [varargout] = likLogistic(hyp, y, mu, s2, inf, i)

% likLogistic - logistic function for binary classification or logit regression.
% The expression for the likelihood is 
%   likLogistic(t) = 1./(1+exp(-t)).
%
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likFunctions.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme. The moments
% \int f^k likLogistic(y,f) N(f|mu,var) df are calculated via a cumulative 
% Gaussian scale mixture approximation.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2013-09-02.
%
% See also likFunctions.m.

if nargin<3, varargout = {'0'}; return; end   % report number of hyperparameters
if nargin>1, y = sign(y); y(y==0) = 1; else y = 1; end % allow only +/- 1 values
if numel(y)==0, y = 1; end

if nargin<5                              % prediction mode if inf is not present
  y = y.*ones(size(mu));                                       % make y a vector
  s2zero = 1; if nargin>3&&numel(s2)>0&&norm(s2)>eps, s2zero = 0; end  % s2==0 ?
  if s2zero                                         % log probability evaluation
    yf = y.*mu;      % product latents and labels
    lp = yf; ok = -35<yf; lp(ok) = -log(1+exp(-yf(ok)));     % log of likelihood
  else                                                              % prediction
    lp = likLogistic(hyp, y, mu, s2, 'infEP');
  end  
  ymu = {}; ys2 = {};
  if nargout>1
    p = exp(lp);
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
      f = mu; yf = y.*f; s = -yf;                   % product latents and labels
      dlp = {}; d2lp = {}; d3lp = {};                         % return arguments
      ps   = max(0,s); 
      lp = -(ps+log(exp(-ps)+exp(s-ps)));                % lp = -(log(1+exp(s)))
      if nargout>1                                           % first derivatives
        s   = min(0,f); 
        p   = exp(s)./(exp(s)+exp(s-f));                    % p = 1./(1+exp(-f))
        dlp = (y+1)/2-p;                          % derivative of log likelihood
        if nargout>2                          % 2nd derivative of log likelihood
          d2lp = -exp(2*s-f)./(exp(s)+exp(s-f)).^2;
          if nargout>3                        % 3rd derivative of log likelihood
            d3lp = 2*d2lp.*(0.5-p);
          end
        end
      end
      varargout = {lp,dlp,d2lp,d3lp};
    else                                                       % derivative mode
      varargout = {[],[],[]};                         % derivative w.r.t. hypers
    end
    
  case 'infEP'
    if nargin<6                                             % no derivative mode
      y = y.*ones(size(mu));                                   % make y a vector
      % likLogistic(t) \approx 1/2 + \sum_{i=1}^5 (c_i/2) erf(lam_i/sqrt(2)t)
      lam = sqrt(2)*[0.44 0.41 0.40 0.39 0.36];    % approx coeffs lam_i and c_i
      c = [1.146480988574439e+02; -1.508871030070582e+03; 2.676085036831241e+03;
          -1.356294962039222e+03;  7.543285642111850e+01                      ];
      [lZc,dlZc,d2lZc] = likErf([], y*ones(1,5), mu*lam, s2*(lam.^2), inf);
      lZ = log_expA_x(lZc,c);       % A=lZc, B=dlZc, d=c.*lam', lZ=log(exp(A)*c)
      dlZ  = expABz_expAx(lZc, c, dlZc, c.*lam');  % ((exp(A).*B)*d)./(exp(A)*c)
      % d2lZ = ((exp(A).*Z)*e)./(exp(A)*c) - dlZ.^2 where e = c.*(lam.^2)'
      d2lZ = expABz_expAx(lZc, c, dlZc.^2+d2lZc, c.*(lam.^2)') - dlZ.^2;
      % The scale mixture approximation does not capture the correct asymptotic
      % behavior; we have linear decay instead of quadratic decay as suggested
      % by the scale mixture approximation. By observing that for large values 
      % of -f*y ln(p(y|f)) for likLogistic is linear in f with slope y, we are
      % able to analytically integrate the tail region.
      val = abs(mu)-196/200*s2-4;       % empirically determined bound at val==0
      lam = 1./(1+exp(-10*val));                         % interpolation weights
      lZtail = min(s2/2-abs(mu),-0.1);  % apply the same to p(y|f) = 1 - p(-y|f)
      dlZtail = -sign(mu); d2lZtail = zeros(size(mu));
      id = y.*mu>0; lZtail(id) = log(1-exp(lZtail(id)));  % label and mean agree
      dlZtail(id) = 0;
      lZ   = (1-lam).*  lZ + lam.*  lZtail;      % interpolate between scale ..
      dlZ  = (1-lam).* dlZ + lam.* dlZtail;              % ..  mixture and   ..
      d2lZ = (1-lam).*d2lZ + lam.*d2lZtail;              % .. tail approximation
      varargout = {lZ,dlZ,d2lZ};
    else                                                       % derivative mode
      varargout = {[]};                                     % deriv. wrt hyp.lik
    end
    
  case 'infVB'
    % variational lower site bound
    % using -log(1+exp(-s)) = s/2 -log( 2*cosh(s/2) );
    % the bound has the form: (b+z/ga)*f - f.^2/(2*ga) - h(ga)/2
    n = numel(s2); b = (y/2).*ones(n,1); z = zeros(size(b));
    varargout = {b,z};
  end
end

%  computes y = log( exp(A)*x ) in a numerically safe way by subtracting the
%  maximal value in each row to avoid cancelation after taking the exp
function y = log_expA_x(A,x)
  N = size(A,2);  maxA = max(A,[],2);      % number of columns, max over columns
  y = log(exp(A-maxA*ones(1,N))*x) + maxA;  % exp(A) = exp(A-max(A))*exp(max(A))
  
%  computes y = ( (exp(A).*B)*z ) ./ ( exp(A)*x ) in a numerically safe way
%  The function is not general in the sense that it yields correct values for
%  all types of inputs. We assume that the values are close together.
function y = expABz_expAx(A,x,B,z)
  N = size(A,2);  maxA = max(A,[],2);      % number of columns, max over columns
  A = A-maxA*ones(1,N);                                 % subtract maximum value
  y = ( (exp(A).*B)*z ) ./ ( exp(A)*x );
