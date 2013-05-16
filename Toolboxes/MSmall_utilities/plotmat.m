function plotmat(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16)

% plotmat( x1, x2, x3,... )
%  or
% plotmat( x1, s1,... )
%
% If size(xi) == [n,2], n>=2, then
% plot( xi(:,1), xi(:,2),... );
%  or
% plot( xi(:,1), xi(:,2), si,... ) if following arg xi is a string
%
% if size(xi) == [2,n] n>2, then 
% plot( xi(1,:), xi(2,:),... );
% etc..
% 

n= 1;

list= [];

while n <= nargin
  
  label= ['x' num2str(n)];
  
  eval(['[j,k]=size(' label ');']);
  
  if j==2 & k~=2
    list= [ list label '(:,1),' label '(:,2)' ];
  else
    list= [ list label '(1,:),' label '(2,:)' ];
  end
  
  if n < nargin
    eval(['s= x' num2str(n+1) ';']);
  else
    s= [];
  end
  if isstr(s)
    n= n+1;
    list= [ list ',x' num2str(n) ];
  end
  
  n= n+1;
  if n<=nargin
    list= [ list ',' ];
  end
  
end

if isstr(list)
  eval( [ 'plot(' list ');' ] );
end