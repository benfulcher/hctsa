function [estimate,nbias,sigma,descriptor] = RM_entropy(x,descriptor,approach,base)
%ENTROPY   Estimates the entropy of stationary signals with
%          independent samples using various approaches.
%   [ESTIMATE,NBIAS,SIGMA,DESCRIPTOR] = ENTROPY(X) or
%   [ESTIMATE,NBIAS,SIGMA,DESCRIPTOR] = ENTROPY(X,DESCRIPTOR) or
%   [ESTIMATE,NBIAS,SIGMA,DESCRIPTOR] = ENTROPY(X,DESCRIPTOR,APPROACH) or
%   [ESTIMATE,NBIAS,SIGMA,DESCRIPTOR] = ENTROPY(X,DESCRIPTOR,APPROACH,BASE)
%
%   ESTIMATE    : The entropy estimate
%   NBIAS       : The N-bias of the estimate
%   SIGMA       : The standard error of the estimate
%   DESCRIPTOR  : The descriptor of the RM_histogram, seel alse ENTROPY
%
%   X           : The time series to be analyzed, a row vector
%   DESCRIPTOR  : Where DESCRIPTOR=[LOWERBOUND,UPPERBOUND,NCELL]
%     LOWERBOUND: Lowerbound of the RM_histogram
%     UPPERBOUND: Upperbound of the RM_histogram
%     NCELL     : The number of cells of the RM_histogram       
%   APPROACH    : The method used, one of the following ones:
%     'unbiased': The unbiased estimate (default)
%     'mmse'    : The minimum mean square error estimate
%     'biased'  : The biased estimate
%   BASE        : The base of the logarithm; default e
%
%   See also: http://www.cs.rug.nl/~rudy/matlab/

%   R. Moddemeijer 
%   Copyright (c) by R. Moddemeijer
%   $Revision: 1.1 $  $Date: 2001/02/05 08:59:36 $


if nargin < 1
   disp('Usage: [ESTIMATE,NBIAS,SIGMA,DESCRIPTOR] = ENTROPY(X)')
   disp('       [ESTIMATE,NBIAS,SIGMA,DESCRIPTOR] = ENTROPY(X,DESCRIPTOR)')
   disp('       [ESTIMATE,NBIAS,SIGMA,DESCRIPTOR] = ENTROPY(X,DESCRIPTOR,APPROACH)')
   disp('       [ESTIMATE,NBIAS,SIGMA,DESCRIPTOR] = ENTROPY(X,DESCRIPTOR,APPROACH,BASE)')
   disp('Where: DESCRIPTOR = [LOWERBOUND,UPPERBOUND,NCELL]')
   return
end

% Some initial tests on the input arguments

[NRowX, NColX] = size(x);

if NRowX~=1
    x = x';
%   error('Invalid dimension of X');
end

if nargin > 4
  error('Too many arguments');
end

if nargin == 1
  [h,descriptor]=RM_histogram(x);
end

if nargin >= 2
  [h,descriptor] = RM_histogram(x,descriptor);
end;

if nargin < 3
  approach = 'unbiased';
end;

if nargin < 4
  base = 2;%exp(1);
end;

lowerbound=descriptor(1);
upperbound=descriptor(2);
ncell=descriptor(3);

estimate=0;
sigma=0;
count=0;

for n = 1:ncell
  if h(n)~=0 
    logf = log(h(n));
  else
    logf = 0;
  end
  count = count+h(n);
  estimate = estimate-h(n)*logf;
  sigma = sigma+h(n)*logf^2;
end

% biased estimate

estimate=estimate/count;
sigma   =sqrt( (sigma/count-estimate^2)/(count-1) );
estimate=estimate+log(count)+log((upperbound-lowerbound)/ncell);
nbias   =-(ncell-1)/(2*count);

% conversion to unbiased estimate

if approach(1)=='u'
  estimate=estimate-nbias;
  nbias=0;
end;

% conversion to minimum mse estimate

if approach(1) == 'm'
  estimate = estimate-nbias;
  nbias = 0;
  lambda = estimate^2/(estimate^2+sigma^2);
  nbias   =(1-lambda)*estimate;
  estimate=lambda*estimate;
  sigma   =lambda*sigma;
end;

% base transformation

estimate=estimate/log(base);
nbias   =nbias   /log(base);
sigma   =sigma   /log(base);

end