% Compute derivative k'(x^p,x^q) of a stationary covariance k(d2) (ard or iso)
% w.r.t. to squared distance d2 = (x^p - x^q)'*inv(P)*(x^p - x^q) measure. Here
% P is either diagonal with ARD parameters ell_1^2,...,ell_D^2 where D is the
% dimension of the input space or ell^2 times the unit matrix for isotropic
% covariance.
% The derivatives can only be computed for one of the following eight
% covariance functions: cov{Matern|PP|RQ|SE}{iso|ard}.
%
% Copyright (c) by Hannes Nickisch, 2013-10-28.
%
% See also INFFITC.M, COVFITC.M.

function Kp = cov_deriv_sq_dist(cov,hyp,x,z)
  if nargin<4, z = []; end                                 % make sure, z exists
  xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;      % determine mode

  if iscell(cov), covstr = cov{1}; else covstr = cov; end
  if ~ischar(covstr), covstr = func2str(covstr); end
  if numel([strfind(covstr,'iso'),strfind(covstr,'ard')])==0
    error('Only iso|ard covariances allowed for derivatives w.r.t. xu.')
  elseif numel([strfind(covstr,'covLIN');
                strfind(covstr,'covGabor');
                strfind(covstr,'covPER')])>0
    error('Gabor|LIN|PER covariances not allowed for derivatives w.r.t. xu.')
  end

  [n,D] = size(x);
  if numel(strfind(covstr,'iso')), id = 1:D; else id = 1; end  % *iso covariance
  ell1 = exp(hyp(1));                        % first characteristic length scale

  Kp = feval(cov{:},hyp,x,z,1);              % use derivative w.r.t. log(ell(1))
  % precompute squared distances
  if dg                                                             % vector kxx
    d2 = zeros(n,1);
  else
    if xeqz                                               % symmetric matrix Kxx
      d2 = sq_dist(x(:,id)'/ell1);
    else                                                 % cross covariances Kxz
      d2 = sq_dist(x(:,id)'/ell1,z(:,id)'/ell1);
    end
  end
  Kp = -1/2*Kp./d2; Kp(d2==0) = 0;