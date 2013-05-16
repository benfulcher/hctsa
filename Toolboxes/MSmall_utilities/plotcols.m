function plotcols(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20)

% plotcols(x1,x2,...) or plotcols(x1,s1,...)
%
% plot the columns of a matrix
% calls plot if dim(xi)==2 or plot3 with first three rows if dim(xi)>2
% either condition must be true for all i
% strings si specify graph parameters


dim= 0;

n= 1;

list= [];

while n <= nargin
  
  label= ['x' num2str(n)];
  
  eval(['[r,c]=size(' label ');']);
  
  if dim~=0 & dim~=r
    error('Not all data of same dimension');
  end
  dim= r;
  
  if dim==1
    list= [ list '1:length(' label '),' label ];
  elseif dim==2
    list= [ list label '(1,:),' label '(2,:)' ];
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
  if dim<=2
    eval( [ 'plot(' list ');' ] );
  else
    eval( [ 'plot3(' list ');' ] );
    xlabel('x1');
    ylabel('x2');
    zlabel('x3');
  end
end
