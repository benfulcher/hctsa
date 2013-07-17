function t=tquant(n, p)
%TQUANT  Quantiles of Student's t distribution
%
%  TQUANT(n, p) is the p-quantile of a t distributed random variable
%  with n degrees of freedom; that is, TQUANT(n, p) is the value below
%  which 100p percent of the t distribution with n degrees of freedom
%  lies.
  
%  Modified 13-Jul-06
%  Author: Tapio Schneider
%          tapio@gps.caltech.edu
%  Nir Krakauer, 19-Mar-11: modified call to fzero to work in Octave

%  References:  
%  L. Devroye, 1986: "Non-Uniform Random Variate Generation", Springer
%  
%  M. Abramowitz and I. A. Stegun, 1964: "Handbook of Mathematical
%     Functions" 
%  
%  See also: tcdf.m in the Matlab Statistics Toolbox (evaluates
%     cumulative distribution function of Student's t)

  if (n ~= round(n) | n < 1)
    error('Usage: TQUANT(n,p) - Degrees of freedom n must be positive integer.')
  end

  if (p<0 | p>1)
    error('Usage: TQUANT(n,p) - Probability p must be in [0,1].')  
  elseif p == 1
    t   = Inf;
    return
  elseif p == 0
    t   = -Inf;
    return
  end
  
  if n == 1
    % Cauchy distribution (cf. Devroye [1986, pp. 29 and 450])
    t   = tan(pi*(p-.5));
  elseif p >= 0.5 
    % positive t-values (cf. M. Abramowitz and I. A. Stegun [1964,
    % Chapter 26])
    b0  = [0, 1];
    f   = @(b) 1 - betainc(b, n/2, .5)/2 - p; 
    b   = fzero(f, b0);    
    t   = sqrt(n/b-n);
  else
    % negative t-values
    b0  = [0, 1];
    f   = @(b) betainc(b, n/2, .5)/2 - p; 
    b   = fzero(f, b0); 
    t   = -sqrt(n/b-n);
  end
