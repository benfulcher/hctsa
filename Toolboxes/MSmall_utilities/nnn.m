function z=nnn(x,n);

% function z=nnn(x,n);
%
% calculates the nth nearest neighbours of the columns of x.
% default n is the sqrt of the number of columns of x.
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk

if nargin < 2
  n=floor(sqrt(length(x(1,:))));
end;

[cx,nx]=size(x);
u= ones(1,nx);

for i=1:nx  
  d = x(:,i)*ones(1,nx) - x;
  if cx==1
    dd= abs(d)';
  else
    dd= sqrt(sum(d.*d))';
  end
  [sorted,ind]=sort(dd);
  z(:,i)= x(:,ind(n));
end;
