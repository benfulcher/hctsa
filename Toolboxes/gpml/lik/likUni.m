function [varargout] = likUni(hyp, y, mu, s2, inf, i)

% likUni - Uniform likelihood function for classification. The expression for 
% the likelihood is 
%   likUni(t) = 1/2.
%
% There are no hyperparameters:
%
% hyp = [ ]
%
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likFunctions.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme.
%
% Copyright (c) by Hannes Nickisch, 2013-09-02.
%
% See also likFunctions.m.

if nargin<3, varargout = {'0'}; return; end   % report number of hyperparameters

if nargin<5                              % prediction mode if inf is not present
  lp = -log(2)*ones(size(mu));
  ymu = {}; ys2 = {};
  if nargout>1
    p = exp(lp);
    ymu = 2*p-1;                                                % first y moment
    if nargout>2
      ys2 = 4*p.*(1-p);                                        % second y moment
    end
  end
  varargout = {lp,ymu,ys2};
else
  switch inf 
  case 'infLaplace'
    if nargin<6                                             % no derivative mode
      lp = -log(2)*ones(size(mu)); dlp = {}; d2lp = {}; d3lp = {};
      if nargout>1
        dlp = zeros(size(mu));               % dlp, derivative of log likelihood
        if nargout>2                    % d2lp, 2nd derivative of log likelihood
          d2lp = zeros(size(mu));
          if nargout>3                  % d3lp, 3rd derivative of log likelihood
            d3lp = zeros(size(mu));
          end
        end
      end
      varargout = {lp,dlp,d2lp,d3lp};
    else                                                       % derivative mode
      varargout = {[],[],[]};                         % derivative w.r.t. hypers
    end

  case 'infEP'
    if nargin<6                                             % no derivative mode
      lZ = -log(2)*ones(size(mu));                           % log part function
      dlZ = {}; d2lZ = {};
      if nargout>1
        dlZ  = zeros(size(mu));                     % 1st derivative w.r.t. mean
        if nargout>2
          d2lZ = zeros(size(mu));                   % 2nd derivative w.r.t. mean
        end
      end
      varargout = {lZ,dlZ,d2lZ};
    else                                                       % derivative mode
      varargout = {[]};                                     % deriv. wrt hyp.lik
    end

  case 'infVB'
    % variational lower site bound
    % t(s) = 1/2
    % the bound has the form: (b+z/ga)*f - f.^2/(2*ga) - h(ga)/2
    n = numel(s2); b = zeros(n,1); y = y.*ones(n,1); z = y;
    varargout = {b,z};
  end
end