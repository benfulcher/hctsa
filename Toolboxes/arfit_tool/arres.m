function [siglev,res, lmp, dof_lmp]=arres(w,A,v,k)
%ARRES 	Test of residuals of fitted AR model.	
%
%  [siglev,res]=ARRES(w,A,v) computes the time series of residuals
% 
%        res(k,:)' = v(k+p,:)'- w - A1*v(k+p-1,:)' - ... - Ap*v(k,:)'
%
%  of an AR(p) model with A=[A1 ... Ap]. If v has three dimensions,
%  the 3rd dimension is treated as indicating multiple realizations
%  (trials) of the time series, that is, v(:,:,itr) is taken as the
%  itr-th trial. 
%       
%  Also returned is the significance level siglev of the modified
%  Li-McLeod portmanteau (LMP) statistic, aggregated over multiple
%  realizations when they are present.
%
%  Correlation matrices for the LMP statistic are computed up to lag
%  k=20, which can be changed to lag k by using
%  [siglev,res]=ARRES(w,A,v,k).

%  Modified 01-Dec-10
%           24-Oct-10 Tim Mullen (added support for multiple realizations)  
%
%  Author: Tapio Schneider
%          tapio@gps.caltech.edu
%
%  Reference:
%    Li, W. K., and A. I. McLeod, 1981: Distribution of the
%        Residual Autocorrelations in Multivariate ARMA Time Series
%        Models, J. Roy. Stat. Soc. B, 43, 231--239.  
%
%  The use of the LMP statistic with correlation matrices of the
%  residuals averaged over the realizations, as we do here, is a
%  heuristic. It may not work well for ntr large. An alternative
%  would be to calculate the LMP statistic for each realization
%  separately and test whether they jointly are consistent with no
%  autocorrelation of the residuals.
  
  n     = size(v,1);                    % number of time steps (per realization)
  m     = size(v,2);                    % dimension of state vectors
  ntr   = size(v,3);                    % number of realizations (trials)
  p     = size(A,2)/m;                  % order of model
  nres  = n-p;                          % number of residuals (per realization)

  % Default value for k 
  if (nargin < 4) 
    k   = min(20, nres-1);
  end
  if (k <= p)                           % check if k is in valid range 
    error('Maximum lag of residual correlation matrices too small.'); 
  end
  if (k >= nres) 
    error('Maximum lag of residual correlation matrices too large.'); 
  end

  w     = w(:)';                        % force w to be row vector

  % Get time series of residuals 
  res = zeros(nres,m,ntr);
  l = 1:nres;                           % vectorized loop l=1,...,nres  
  res(l,:,:) = v(l+p,:,:) - repmat(w, [nres, 1, ntr]);
  for itr=1:ntr
    for j=1:p
      res(l,:,itr) = res(l,:,itr) - v(l-j+p,:,itr)*A(:, (j-1)*m+1:j*m)';
    end
  end
  % end of loop over l
  
  % For computation of correlation matrices, center residuals by
  % subtraction of the mean  
  resc  = res - repmat(mean(res), [nres, 1, 1]);
  
  % Compute correlation matrix of the residuals
  % compute lag zero correlation matrix
  c0    = zeros(m, m);
  for itr=1:ntr
    c0  = c0 + resc(:,:,itr)'*resc(:,:,itr);
  end
  d     = diag(c0);
  dd    = sqrt(d*d');
  c0    = c0./dd;

  % compute lag l correlation matrix
  cl    = zeros(m, m, k);
  for l=1:k
    for itr=1:ntr
      cl(:,:,l) = cl(:,:,l) ...
          + resc(1:nres-l, :, itr)'*resc(l+1:nres, :, itr);
    end
    cl(:,:,l) = cl(:,:,l)./dd;  
  end
  
  % Get "covariance matrix" in LMP statistic
  c0_inv= inv(c0);                      % inverse of lag 0 correlation matrix
  rr    = kron(c0_inv, c0_inv);         % "covariance matrix" in LMP statistic

  % Compute modified Li-McLeod portmanteau statistic
  lmp   = 0;                            % LMP statistic initialization
  x     = zeros(m*m,1);                 % correlation matrix arranged as vector
  for l=1:k
    x   = reshape(cl(:,:,l), m^2, 1);   % arrange cl as vector by stacking columns
    lmp = lmp + x'*rr*x;                % sum up LMP statistic
  end
  ntot  = nres*ntr;                     % total number of residual vectors
  lmp   = ntot*lmp + m^2*k*(k+1)/2/ntot;% add remaining term and scale
  dof_lmp = m^2*(k-p);                  % degrees of freedom for LMP statistic
      
  % Significance level with which hypothesis of uncorrelatedness is rejected
  siglev = 1 - gammainc(lmp/2, dof_lmp/2);