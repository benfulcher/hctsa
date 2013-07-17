function [varargout] = likGauss(hyp, y, mu, s2, inf, i)

% likGauss - Gaussian likelihood function for regression. The expression for the 
% likelihood is 
%   likGauss(t) = exp(-(t-y)^2/2*sn^2) / sqrt(2*pi*sn^2),
% where y is the mean and sn is the standard deviation.
%
% The hyperparameters are:
%
% hyp = [  log(sn)  ]
%
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likFunctions.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2013-01-21
%
% See also LIKFUNCTIONS.M.

if nargin<3, varargout = {'1'}; return; end   % report number of hyperparameters

sn2 = exp(2*hyp);

if nargin<5                              % prediction mode if inf is not present
  if numel(y)==0,  y = zeros(size(mu)); end
  s2zero = 1; if nargin>3, if norm(s2)>0, s2zero = 0; end, end         % s2==0 ?
  if s2zero                                                    % log probability
    lp = -(y-mu).^2./sn2/2-log(2*pi*sn2)/2; s2 = 0;
  else
    lp = likGauss(hyp, y, mu, s2, 'infEP');                         % prediction
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
    if nargin<6                                             % no derivative mode
      if numel(y)==0, y=0; end
      ymmu = y-mu; dlp = {}; d2lp = {}; d3lp = {};
      lp = -ymmu.^2/(2*sn2) - log(2*pi*sn2)/2; 
      if nargout>1
        dlp = ymmu/sn2;                      % dlp, derivative of log likelihood
        if nargout>2                    % d2lp, 2nd derivative of log likelihood
          d2lp = -ones(size(ymmu))/sn2;
          if nargout>3                  % d3lp, 3rd derivative of log likelihood
            d3lp = zeros(size(ymmu));
          end
        end
      end
      varargout = {lp,dlp,d2lp,d3lp};
    else                                                       % derivative mode
      lp_dhyp = (y-mu).^2/sn2 - 1;  % derivative of log likelihood w.r.t. hypers
      dlp_dhyp = 2*(mu-y)/sn2;                               % first derivative,
      d2lp_dhyp = 2*ones(size(mu))/sn2;   % and also of the second mu derivative
      varargout = {lp_dhyp,dlp_dhyp,d2lp_dhyp};
    end

  case 'infEP'
    if nargin<6                                             % no derivative mode
      lZ = -(y-mu).^2./(sn2+s2)/2 - log(2*pi*(sn2+s2))/2;    % log part function
      dlZ = {}; d2lZ = {};
      if nargout>1
        dlZ  = (y-mu)./(sn2+s2);                    % 1st derivative w.r.t. mean
        if nargout>2
          d2lZ = -1./(sn2+s2);                      % 2nd derivative w.r.t. mean
        end
      end
      varargout = {lZ,dlZ,d2lZ};
    else                                                       % derivative mode
      dlZhyp = ((y-mu).^2./(sn2+s2)-1) ./ (1+s2./sn2);   % deriv. w.r.t. hyp.lik
      varargout = {dlZhyp};
    end

  case 'infVB'
    if nargin<6
      % variational lower site bound
      % t(s) = exp(-(y-s)^2/2sn2)/sqrt(2*pi*sn2)
      % the bound has the form: b*s - s.^2/(2*ga) - h(ga)/2 with b=y/ga
      ga = s2; n = numel(ga); b = y./ga; y = y.*ones(n,1);
      db = -y./ga.^2; d2b = 2*y./ga.^3;
      h = zeros(n,1); dh = h; d2h = h;         % allocate memory for return args
      id = ga(:)<=sn2+1e-8;                            % OK below noise variance
      h(id) = y(id).^2./ga(id) + log(2*pi*sn2); h(~id) = Inf;
      dh(id) = -y(id).^2./ga(id).^2;
      d2h(id) = 2*y(id).^2./ga(id).^3;
      id = ga<0; h(id) = Inf; dh(id) = 0; d2h(id) = 0;     % neg. var. treatment
      varargout = {h,b,dh,db,d2h,d2b};
    else
      ga = s2; n = numel(ga); 
      dhhyp = zeros(n,1); dhhyp(ga(:)<=sn2) = 2;
      dhhyp(ga<0) = 0;              % negative variances get a special treatment
      varargout = {dhhyp};                               % deriv. w.r.t. hyp.lik
    end
  end
end
