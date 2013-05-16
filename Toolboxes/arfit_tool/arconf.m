function [Aerr, werr]=arconf(A, C, w, th)
%ARCONF	Confidence intervals for AR coefficients.
%
%  For an AR(p) model that has been fitted with ARFIT,
%  [Aerr,werr]=ARCONF(A,C,w,th) computes the margins of error Aerr and
%  werr such that (A +/- Aerr) and (w +/- werr) are approximate 95%
%  confidence intervals for the elements of the coefficient matrix A
%  and for the components of the intercept vector w. The input
%  arguments of ARCONF are output of AR.
%
%  If no intercept vector w has been fitted with ARFIT (i.e., the flag
%  'zero' was an input argument of ARFIT), then [Aerr]=ARCONF(A,C,th)
%  computes the margins of error only for the elements of the
%  coefficient matrix A.
%
%  The confidence intervals are based on Student's t distribution,
%  which for small samples yields only approximate confidence
%  intervals. Inferences drawn from small samples must therefore be
%  interpreted cautiously.
%
%  See also ARFIT.

%  Modified 30-Dec-99
%  Author: Tapio Schneider
%          tapio@gps.caltech.edu

  ccoeff = .95;            % confidence coefficient
  m 	 = size(C,1);      % dimension of state space
  p 	 = size(A,2)/m;    % order of model

  if (nargin == 3)
    %  no intercept vector has been fitted
    Aaug = A;
    th 	 = w;
    w 	 = [];
    np 	 = m*p;            % number of parameter vectors of size m
  else
    Aaug = [w A];
    np 	 = m*p+1;          % number of parameter vectors of size m
  end
  % number of degrees of freedom for residual covariance matrix
  dof 	 = th(1,1);               
  % quantile of t distribution for given confidence coefficient and dof
  t      = tquant(dof, .5+ccoeff/2);
  
  % Get matrix Uinv that appears in the covariance matrix of the least squares
  % estimator
  Uinv   = th(2:size(th,1), :);

  % Compute approximate confidence intervals for elements of Aaug  
  Aaug_err = zeros(m, np);
  for j=1:m
    for k=1:np
      Aaug_err(j,k) = t * sqrt( Uinv(k ,k)* C(j,j) );
    end
  end

  if (nargin == 3)
    %  No intercept vector has been fitted
    Aerr  = Aaug_err;
  else 
    % An intercept vector has been fitted => return margins of error
    % for intercept vector and for AR coefficients separately
    werr  = Aaug_err(:, 1);
    Aerr  = Aaug_err(:, 2:np);
  end
