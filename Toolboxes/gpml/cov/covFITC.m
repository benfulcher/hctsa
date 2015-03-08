function [K,Kuu,Ku] = covFITC(cov, xu, hyp, x, z, i)

% Covariance function to be used together with the FITC approximation.
%
% The function allows for more than one output argument and does not respect the
% interface of a proper covariance function. In fact, it wraps a proper
% covariance function such that it can be used together with infFITC.m.
% Instead of outputing the full covariance, it returns cross-covariances between
% the inputs x, z and the inducing inputs xu as needed by infFITC.m
%
% Copyright (c) by Ed Snelson, Carl Edward Rasmussen 
%                                               and Hannes Nickisch, 2011-11-02.
%
% See also COVFUNCTIONS.M, INFFITC.M.

if nargin<4, K = feval(cov{:}); return, end
if nargin<5, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

if size(xu,2) ~= size(x,2)
  error('Dimensionality of inducing inputs must match training inputs');
end

if nargin<6                                                        % covariances
  if dg
    K = feval(cov{:},hyp,x,'diag');
  else
    if xeqz
        K = feval(cov{:},hyp,x,'diag');
        if nargout>1, Kuu = feval(cov{:},hyp,xu);   end
        if nargout>2, Ku  = feval(cov{:},hyp,xu,x); end
    else
      K = feval(cov{:},hyp,xu,z);
    end
  end
else                                                               % derivatives
  if dg
    K = feval(cov{:},hyp,x,'diag',i);
  else
    if xeqz
        K   = feval(cov{:},hyp,x,'diag',i);
        if nargout>1, Kuu = feval(cov{:},hyp,xu,[],i); end
        if nargout>2, Ku  = feval(cov{:},hyp,xu,x,i);  end
    else
      K = feval(cov{:},hyp,xu,z,i);
    end
  end
end