function [siglev,res]=arres(w,A,v,k)
%ARRES 	Test of residuals of fitted AR model.	
%
%  [siglev,res]=ARRES(w,A,v) computes the time series of residuals
% 
%        res(k,:)' = v(k+p,:)'- w - A1*v(k+p-1,:)' - ... - Ap*v(k,:)'
%
%  of an AR(p) model with A=[A1 ... Ap]. 
%       
%  Also returned is the significance level siglev of the modified
%  Li-McLeod portmanteau (LMP) statistic.
%
%  Correlation matrices for the LMP statistic are computed up to lag
%  k=20, which can be changed to lag k by using
%  [siglev,res]=ARRES(w,A,v,k).

%  Modified 17-Dec-99
%  Author: Tapio Schneider
%          tapio@gps.caltech.edu
%
%  Reference:
%    Li, W. K., and A. I. McLeod, 1981: Distribution of the
%        Residual Autocorrelations in Multivariate ARMA Time Series
%        Models, J. Roy. Stat. Soc. B, 43, 231--239.  
  
  m     = size(v,2);                    % dimension of state vectors
  p     = size(A,2)/m;                  % order of model
  n     = length(v);                    % number of observations
  nres  = n-p;                          % number of residuals

  % Default value for k 
  if (nargin < 4) 
    k   = 20;
  end
  if (k <= p)                           % check if k is in valid range 
    error('Maximum lag of residual correlation matrices too small.'); 
  end
  if (k >= nres) 
    error('Maximum lag of residual correlation matrices too large.'); 
  end

  w     = w(:)';                        % force w to be row vector

  % Get time series of residuals 
  l = 1:nres; 				% vectorized loop l=1,...,nres 
  res(l,:) = v(l+p,:) - ones(nres,1)*w;
  for j=1:p
    res(l,:) = res(l,:) - v(l-j+p,:)*A(:, (j-1)*m+1:j*m)';
  end
  % end of loop over l
  
  % Center residuals by subtraction of the mean 
  res   = res - ones(nres,1)*mean(res);
  
  % Compute lag zero correlation matrix of the residuals
  c0    = res'*res;
  d     = diag(c0);
  dd    = sqrt(d*d');
  c0    = c0./dd;
  
  % Get "covariance matrix" in LMP statistic
  c0_inv= inv(c0);                      % inverse of lag 0 correlation matrix
  rr    = kron(c0_inv, c0_inv);         % "covariance matrix" in LMP statistic

  % Initialize LMP statistic and correlation matrices
  lmp   = 0;                            % LMP statistic
  cl    = zeros(m,m);                   % correlation matrix
  x     = zeros(m*m,1);                 % correlation matrix arranged as vector
  
  % Compute modified Li-McLeod portmanteau statistic
  for l=1:k
    cl  = (res(1:nres-l, :)'*res(l+1:nres,:))./dd;  % lag l correlation matrix
    x   = reshape(cl,m*m, 1);           % arrange cl as vector by stacking columns
    lmp = lmp + x'*rr*x;                % sum up LMP statistic
  end
  lmp   = n*lmp + m^2*k*(k+1)/2/n;      % add remaining term and scale
  dof_lmp = m^2*(k-p);                  % degrees of freedom for LMP statistic

  % Significance level with which hypothesis of uncorrelatedness is rejected
  siglev = 1 - gammainc(lmp/2, dof_lmp/2);







