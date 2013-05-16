function [w,s] = normalize(v)

% [w,s] = normalize(v)
%
% normalize the columns of v
% s is the row vector of normalizing factors
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk

[n,m]= size(v);
w= zeros(n,m);
s= ones(1,m);
for i=1:m
  n= norm(v(:,i));
  if n~=0
    w(:,i)= v(:,i) / n;
    s(i)= n;
  end
end
