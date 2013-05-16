function p=plotmat3(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20)

% plotmat3( x1, x2, x3,... )
%  or
% plotmat3( x1, s1,... )
%
% If size(xi) == [n,3], then
% plot3( xi(:,1), xi(:,2),xi(:,3),... );
%  or
% plot3( xi(:,1), xi(:,2),xi(:,3), si,... ) if following arg xi is a string
%
% if size(xi) == [3,n] n~=3, then 
% plot3( xi(1,:), xi(2,:),xi(3,:),... );
% etc..
% 

n= 1;

list= [];

while n <= nargin
  
  label= ['x' num2str(n)];
  
  eval(['[j,k]=size(' label ');']);
  
  if j==3 & k~=3
    list= [ list label '(:,1),' label '(:,2),' label '(:,3)' ];
  else
    list= [ list label '(1,:),' label '(2,:),' label '(3,:)' ];
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
  eval( [ 'p=plot3(' list ');' ] );
end
