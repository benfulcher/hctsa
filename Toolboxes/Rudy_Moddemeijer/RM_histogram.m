function [result,descriptor] = RM_histogram(x,descriptor)
% HISTOGRAM Computes the frequency histogram of the row vector x.
%   [RESULT,DESCRIPTOR] = HISTOGRAM(X) or
%   [RESULT,DESCRIPTOR] = HISTOGRAM(X,DESCRIPTOR) or
% where
%   DESCRIPTOR = [LOWER,UPPER,NCELL]
%
%   RESULT    : A row vector containing the histogram
%   DESCRIPTOR: The used descriptor
%
%   X         : The row vector be analyzed
%   DESCRIPTOR: The descriptor of the histogram
%     LOWER   : The lowerbound of the histogram
%     UPPER   : The upperbound of the histogram
%     NCELL   : The number of cells of the histogram
%
%   See also: http://www.cs.rug.nl/~rudy/matlab/

%   R. Moddemeijer 
%   Copyright (c) by R. Moddemeijer
%   $Revision: 1.2 $  $Date: 2001/02/05 09:54:29 $


if nargin < 1
   disp('Usage: RESULT = HISTOGRAM(X)')
   disp('       RESULT = HISTOGRAM(X,DESCRIPTOR)')
   disp('Where: DESCRIPTOR = [LOWER,UPPER,NCELL]')
   return
end

% Some initial tests on the input arguments

[NRowX, NColX] = size(x);

if NRowX ~= 1
  error('Invalid dimension of X');
end;

if nargin > 2
  error('Too many arguments');
end;

if nargin==1
  minx=min(x);
  maxx=max(x);
  delta=(maxx-minx)/(length(x)-1);
  ncell=ceil(sqrt(length(x)));
  descriptor=[minx-delta/2,maxx+delta/2,ncell];
end;

lower=descriptor(1);
upper=descriptor(2);
ncell=descriptor(3);

if ncell<1 
  error('Invalid number of cells')
end;

if upper<=lower
  error('Invalid bounds')
end;

result(1:ncell)=0;

y=round( (x-lower)/(upper-lower)*ncell + 1/2 );
for n=1:NColX
  index=y(n);
  if index >= 1 & index<=ncell
    result(index)=result(index)+1;
  end;
end;

