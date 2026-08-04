function [K,dK] = covPref(cov, hyp, x, z)

% covPref - covariance function for preference learning. The covariance
% function corresponds to a prior on f(x1) - f(x2).
%
% k(x,z) = k_0(x1,z1) + k_0(x2,z2) - k_0(x1,z2) - k_0(x2,z1).
%
% The hyperparameters are:
%
% hyp = [ hyp_k0 ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% See Collaborative Gaussian Processes for Preference Learning, NIPS 2014.
%
% Copyright (c) by Hannes Nickisch and Roman Garnett, 2016-04-17.
%
% See also covFunctions.m.

if nargin<3, K = strrep(feval(cov{:}),'D','D/2'); return; end     % no of params
if nargin<4, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

x1 = x(:,1:end/2); x2 = x(:,1+end/2:end);
if dg || xeqz
  z1 = x1; z2 = x2;
else
  z1 = z(:,1:end/2); z2 = z(:,1+end/2:end);
end
if xeqz
  [K11,dK11] = feval(cov{:},hyp,x1);
  [K22,dK22] = feval(cov{:},hyp,x2);
else
  if dg
    [K11,dK11] = feval(cov{:},hyp,x1,'diag');
    [K22,dK22] = feval(cov{:},hyp,x2,'diag');
  else
    [K11,dK11] = feval(cov{:},hyp,x1,z1);
    [K22,dK22] = feval(cov{:},hyp,x2,z2);  
  end
end
[K12,dK12] = feval(cov{:},hyp,x1,z2);
[K21,dK21] = feval(cov{:},hyp,x2,z1);

if dg, K12 = diag(K12); K21 = diag(K21); end
K = K11 + K22 - K12 - K21;
dK = @(Q) dirder(Q,dK11,dK22,dK12,dK21,dg,xeqz);

function [dhyp,dx] = dirder(Q,dK11,dK22,dK12,dK21,dg,xeqz)
  if nargout > 1
    [dhyp11,dx11] = dK11(Q); [dhyp22,dx22] = dK22(Q);
    if dg
      [dhyp12,dx12] = dK12(diag(Q)); [dhyp21,dx21] = dK21(diag(Q));
    else
      [dhyp12,dx12] = dK12(Q); [dhyp21,dx21] = dK21(Q);
    end
    if xeqz
      [junk,dx21t] = dK12(Q'); [junk,dx12t] = dK21(Q');
      dx = [dx11-dx12-dx21t, dx22-dx21-dx12t];
    else
      dx = [dx11-dx12, dx22-dx21];  
    end
  else
    dhyp11 = dK11(Q); dhyp22 = dK22(Q);
    if dg
      dhyp12 = dK12(diag(Q)); dhyp21 = dK21(diag(Q)); 
    else
      dhyp12 = dK12(Q); dhyp21 = dK21(Q);
    end   
  end
  dhyp = dhyp11 + dhyp22 - dhyp12 - dhyp21;