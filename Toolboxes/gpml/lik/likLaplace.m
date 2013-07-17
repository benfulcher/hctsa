function [varargout] = likLaplace(hyp, y, mu, s2, inf, i)

% likLaplace - Laplacian likelihood function for regression. 
% The expression for the likelihood is 
%   likLaplace(t) = exp(-|t-y|/b)/(2*b) with b = sn/sqrt(2),
% where y is the mean and sn^2 is the variance.
%
% The hyperparameters are:
%
% hyp = [  log(sn)  ]
%
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likFunctions.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme. 
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2012-11-04.
%
% See also LIKFUNCTIONS.M.

if nargin<3, varargout = {'1'}; return; end   % report number of hyperparameters

sn = exp(hyp); b = sn/sqrt(2);
if nargin<5                              % prediction mode if inf is not present
  if numel(y)==0, y = zeros(size(mu)); end
  s2zero = 1; if nargin>3, if norm(s2)>0, s2zero = 0; end, end         % s2==0 ?
  if s2zero                                         % log probability evaluation
    lp = -abs(y-mu)./b -log(2*b); s2 = 0;
  else                                                              % prediction
    lp = likLaplace(hyp, y, mu, s2, 'infEP');
  end
  ymu = {}; ys2 = {};
  if nargout>1
    ymu = mu;                                                   % first y moment
    if nargout>2
      ys2 = s2 + sn.^2;                                        % second y moment
    end
  end
  varargout = {lp,ymu,ys2};
else                                                            % inference mode
  switch inf 
  case 'infLaplace'
    if nargin<6                                             % no derivative mode
      if numel(y)==0, y=0; end
      ymmu = y-mu; dlp = {}; d2lp = {}; d3lp = {};     
      lp = -abs(ymmu)/b -log(2*b);
      if nargout>1
        dlp = sign(ymmu)/b;                  % dlp, derivative of log likelihood
        if nargout>2                    % d2lp, 2nd derivative of log likelihood
          d2lp = zeros(size(ymmu));
          if nargout>3                  % d3lp, 3rd derivative of log likelihood
            d3lp = zeros(size(ymmu));
          end
        end
      end
      varargout = {lp,dlp,d2lp,d3lp};
    else                                                       % derivative mode
      lp_dhyp = abs(y-mu)/b - 1;    % derivative of log likelihood w.r.t. hypers
      dlp_dhyp = sign(mu-y)/b;                               % first derivative,
      d2lp_dhyp = zeros(size(mu));        % and also of the second mu derivative
      varargout = {lp_dhyp,dlp_dhyp,d2lp_dhyp};
    end
    
  case 'infEP'
    n = max([numel(y),numel(mu),numel(s2),numel(sn)]); on = ones(n,1);
    y = y(:).*on; mu = mu(:).*on; s2 = s2(:).*on; sn = sn(:).*on; % vectors only
    fac = 1e3;          % factor between the widths of the two distributions ...
       % ... from when one considered a delta peak, we use 3 orders of magnitude
    idlik = fac*sn<sqrt(s2);                        % Likelihood is a delta peak
    idgau = fac*sqrt(s2)<sn;                          % Gaussian is a delta peak
    id = ~idgau & ~idlik;                          % interesting case in between
    if nargin<6                                             % no derivative mode
      lZ = zeros(n,1); dlZ = lZ; d2lZ = lZ;                    % allocate memory
      if any(idlik)
        [lZ(idlik),dlZ(idlik),d2lZ(idlik)] = ...
                                likGauss(log(s2(idlik))/2, mu(idlik), y(idlik));
      end
      if any(idgau)
        [lZ(idgau),dlZ(idgau),d2lZ(idgau)] = ...
                                likLaplace(log(sn(idgau)), mu(idgau), y(idgau));
      end
      if any(id)

        % substitution to obtain unit variance, zero mean Laplacian
        tmu = (mu(id)-y(id))./sn(id); tvar = s2(id)./sn(id).^2;

%         % an equivalent implementation based on lerfc(t) = log(erfc(t)) 
%         % instead of logphi; we found it to be less stable and therefore
%         % use logphi
%         zp  = (tvar+tmu/sqrt(2))./sqrt(tvar);  vp = tvar+sqrt(2)*tmu;
%         zm  = (tvar-tmu/sqrt(2))./sqrt(tvar);  vm = tvar-sqrt(2)*tmu;
%         lM0p = lerfc(zp);                                         % 0th moment
%         lM0m = lerfc(zm);
%         lZ(id) = logsum2exp([vp+lM0p, vm+lM0m]) - log(2*sqrt(2)*sn(id));
%         if nargout>1                                              % 1st moment
%           M1p  = sqrt(tvar).*( zp.*exp(lM0p) - exp(-zp.*zp)/sqrt(pi) );
%           M1m  = sqrt(tvar).*( zm.*exp(lM0m) - exp(-zm.*zm)/sqrt(pi) );
%           lZmax = max([vp-lZ(id),vm-lZ(id)],[],2);
%           m1m0 = exp(lZmax).*( exp(vp-lZ(id)-lZmax).*M1p ...
%                               -exp(vm-lZ(id)-lZmax).*M1m )/2 + y(id);
%           dlZ(id) = (m1m0-mu(id))./s2(id);
%           if nargout>2                                            % 2nd moment
%             M2p = 2*tvar.*exp(vp-lZ(id)+lM0p)/2+(2*tvar+sqrt(2)*tmu) ...
%                                                   .*M1p.*exp(vp-lZ(id));
%             M2m = 2*tvar.*exp(vm-lZ(id)+lM0m)/2+(2*tvar-sqrt(2)*tmu) ...
%                                                   .*M1m.*exp(vm-lZ(id));
%             m2m0 = (M2p+M2m).*sn(id)/sqrt(8) + 2*y(id).*m1m0 - y(id).^2;
%             d2lZ(id) = (m2m0 - m1m0.^2)./s2(id).^2 - 1./s2(id);
%           end
%         end

        % an implementation based on logphi(t) = log(normcdf(t))
        zp = (tmu+sqrt(2)*tvar)./sqrt(tvar);
        zm = (tmu-sqrt(2)*tvar)./sqrt(tvar);
        ap =  logphi(-zp)+sqrt(2)*tmu;
        am =  logphi( zm)-sqrt(2)*tmu;
        lZ(id) = logsum2exp([ap,am]) + tvar - log(sn(id)*sqrt(2));
        
        if nargout>1
          lqp = -zp.^2/2 - log(2*pi)/2 - logphi(-zp);       % log( N(z)/Phi(z) )
          lqm = -zm.^2/2 - log(2*pi)/2 - logphi( zm);
          dap = -exp(lqp-log(s2(id))/2) + sqrt(2)./sn(id);
          dam =  exp(lqm-log(s2(id))/2) - sqrt(2)./sn(id);
                        % ( exp(ap).*dap + exp(am).*dam )./( exp(ap) + exp(am) )
          dlZ(id) = expABz_expAx([ap,am],[1;1],[dap,dam],[1;1]);
          
          if nargout>2
            a = sqrt(8)./sn(id)./sqrt(s2(id));
            bp = 2./sn(id).^2 - (a - zp./s2(id)).*exp(lqp);
            bm = 2./sn(id).^2 - (a + zm./s2(id)).*exp(lqm);
            % d2lZ(id) = ( exp(ap).*bp + exp(am).*bm )./( exp(ap) + exp(am) ) ...
            %            - dlZ(id).^2;
            d2lZ(id) = expABz_expAx([ap,am],[1;1],[bp,bm],[1;1]) - dlZ(id).^2;
          end
        end

      end
      varargout = {lZ,dlZ,d2lZ};
    else                                                       % derivative mode
      dlZhyp = zeros(n,1);
      if any(idlik)
        dlZhyp(idlik) = 0;
      end
      if any(idgau)
        dlZhyp(idgau) = ...
               likLaplace(log(sn(idgau)), mu(idgau), y(idgau), 'infLaplace', 1);
      end
      if any(id)
        % substitution to obtain unit variance, zero mean Laplacian
        tmu = (mu(id)-y(id))./sn(id);        tvar = s2(id)./sn(id).^2;
        zp  = (tvar+tmu/sqrt(2))./sqrt(tvar);  vp = tvar+sqrt(2)*tmu;
        zm  = (tvar-tmu/sqrt(2))./sqrt(tvar);  vm = tvar-sqrt(2)*tmu;
        dzp = (-s2(id)./sn(id)+tmu.*sn(id)/sqrt(2)) ./ sqrt(s2(id));
        dvp = -2*tvar - sqrt(2)*tmu;
        dzm = (-s2(id)./sn(id)-tmu.*sn(id)/sqrt(2)) ./ sqrt(s2(id));
        dvm = -2*tvar + sqrt(2)*tmu;
        lezp = lerfc(zp); % ap = exp(vp).*ezp
        lezm = lerfc(zm); % am = exp(vm).*ezm
        vmax = max([vp+lezp,vm+lezm],[],2); % subtract max to avoid numerical pb
        ep  = exp(vp+lezp-vmax);
        em  = exp(vm+lezm-vmax);
        dap = ep.*(dvp - 2/sqrt(pi)*exp(-zp.^2-lezp).*dzp);
        dam = em.*(dvm - 2/sqrt(pi)*exp(-zm.^2-lezm).*dzm);        
        dlZhyp(id) = (dap+dam)./(ep+em) - 1;       
      end
      varargout = {dlZhyp};                                  % deriv. wrt hypers
    end

  case 'infVB'
    if nargin<6
      % variational lower site bound
      % t(s) = exp(-sqrt(2)|y-s|/sn) / sqrt(2*sn²)
      % the bound has the form: b*s - s.^2/(2*ga) - h(ga)/2 with b=y/ga!!
      ga = s2; n = numel(ga); b = y./ga; y = y.*ones(n,1);
      db = -y./ga.^2; d2b = 2*y./ga.^3;
      h   = 2*ga/sn^2 + log(2*sn^2) + y.^2./ga;
      dh  = 2/sn^2 - y.^2./ga.^2;
      d2h = 2*y.^2./ga.^3;
      id = ga<0; h(id) = Inf; dh(id) = 0; d2h(id) = 0;     % neg. var. treatment
      varargout = {h,b,dh,db,d2h,d2b};
    else
      ga = s2; dhhyp = -4*ga/sn^2 + 2;
      dhhyp(ga<0) = 0;              % negative variances get a special treatment
      varargout = {dhhyp};                                  % deriv. wrt hyp.lik
    end
  end
end

% computes y = log( sum(exp(x),2) ) in a numerically safe way by subtracting 
%  the row maximum to avoid cancelation after taking the exp
%  the sum is done along the rows
function [y,x] = logsum2exp(logx)
  N = size(logx,2);
  max_logx = max(logx,[],2);
  % we have all values in the log domain, and want to calculate a sum
  x = exp(logx-max_logx*ones(1,N));
  y = log(sum(x,2)) + max_logx;

% numerically safe implementation of f(t) = log(1-erf(t)) = log(erfc(t))
function f = lerfc(t)
  f  = zeros(size(t));
  tmin = 20; tmax = 25;
  ok = t<tmin;                               % log(1-erf(t)) is safe to evaluate
  bd = t>tmax;                                            % evaluate tight bound
  interp = ~ok & ~bd;                % linearly interpolate between both of them
  f(~ok) = log(2/sqrt(pi)) -t(~ok).^2 -log(t(~ok)+sqrt( t(~ok).^2+4/pi ));
  lam = 1./(1+exp( 12*(1/2-(t(interp)-tmin)/(tmax-tmin)) ));   % interp. weights
  f(interp) = lam.*f(interp) + (1-lam).*log(erfc( t(interp) ));
  f(ok) = f(ok) + log(erfc( t(ok) ));                                % safe eval
%  computes y = ( (exp(A).*B)*z ) ./ ( exp(A)*x ) in a numerically safe way
%  The function is not general in the sense that it yields correct values for
%  all types of inputs. We assume that the values are close together.

function y = expABz_expAx(A,x,B,z)
  N = size(A,2);  maxA = max(A,[],2);      % number of columns, max over columns
  A = A-maxA*ones(1,N);                                 % subtract maximum value
  y = ( (exp(A).*B)*z ) ./ ( exp(A)*x ); 

% safe implementation of the log of phi(x) = \int_{-\infty}^x N(f|0,1) df
% logphi(z) = log(normcdf(z))
function lp = logphi(z)
  lp = zeros(size(z));                                         % allocate memory
  zmin = -6.2; zmax = -5.5;
  ok = z>zmax;                                % safe evaluation for large values
  bd = z<zmin;                                                 % use asymptotics
  ip = ~ok & ~bd;                             % interpolate between both of them
  lam = 1./(1+exp( 25*(1/2-(z(ip)-zmin)/(zmax-zmin)) ));       % interp. weights
  lp( ok) = log( (1+erf(z(ok)/sqrt(2)))/2 );
  % use lower and upper bound acoording to Abramowitz&Stegun 7.1.13 for z<0
  % lower -log(pi)/2 -z.^2/2 -log( sqrt(z.^2/2+2   ) -z/sqrt(2) )
  % upper -log(pi)/2 -z.^2/2 -log( sqrt(z.^2/2+4/pi) -z/sqrt(2) )
  % the lower bound captures the asymptotics
  lp(~ok) = -log(pi)/2 -z(~ok).^2/2 -log( sqrt(z(~ok).^2/2+2)-z(~ok)/sqrt(2) );
  lp( ip) = (1-lam).*lp(ip) + lam.*log( (1+erf(z(ip)/sqrt(2)))/2 );
