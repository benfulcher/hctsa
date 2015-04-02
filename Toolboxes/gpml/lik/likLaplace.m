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
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2013-10-16.
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

        % an implementation based on logphi(t) = log(normcdf(t))
        zp = (tmu+sqrt(2)*tvar)./sqrt(tvar);
        zm = (tmu-sqrt(2)*tvar)./sqrt(tvar);
        ap =  logphi(-zp)+sqrt(2)*tmu;
        am =  logphi( zm)-sqrt(2)*tmu;
        lZ(id) = logsumexp2([ap,am]) + tvar - log(sn(id)*sqrt(2));
        
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
            % d2lZ(id) = ( exp(ap).*bp + exp(am).*bm )./( exp(ap) + exp(am) )...
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
        lezp = logerfc(zp); % ap = exp(vp).*ezp
        lezm = logerfc(zm); % am = exp(vm).*ezm
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
    % variational lower site bound
    % t(s) = exp(-sqrt(2)|y-s|/sn) / sqrt(2*sn²)
    % the bound has the form: (b+z/ga)*f - f.^2/(2*ga) - h(ga)/2
    n = numel(s2); b = zeros(n,1); y = y.*ones(n,1); z = y;
    varargout = {b,z};
  end
end

% logerfc(z) = log(1-erf(z))
function lc = logerfc(z)
  lc = logphi(-z*sqrt(2)) + log(2);

function y = expABz_expAx(A,x,B,z)
  N = size(A,2);  maxA = max(A,[],2);      % number of columns, max over columns
  A = A-maxA*ones(1,N);                                 % subtract maximum value
  y = ( (exp(A).*B)*z ) ./ ( exp(A)*x );