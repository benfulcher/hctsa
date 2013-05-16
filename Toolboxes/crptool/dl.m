function [a_out, b_out]=dl(x)
% DL   Mean of the diagonal line lengths and their distribution.
%    A=DL(X) computes the mean of the length of the diagonal 
%    line structures in a recurrence plot.
%
%    [A B]=DL(X) computes the mean A and the lengths of the
%    found diagonal lines, stored in B. In order to get the 
%    histogramme of the line lengths, simply call 
%    HIST(B,[1 MAX(B)]).
%
%    Examples: X = crp(rand(200,1),1,1,.3,'fan','silent');
%              [l l_dist] = dl(X);
%              hist(l_dist,200)
%
%    See also CRQA, TT.

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 2001-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:32:09 $
% $Revision: 3.4 $
%
% $Log: dl.m,v $
% Revision 3.4  2009/03/24 08:32:09  marwan
% copyright address changed
%
% Revision 3.3  2005/11/23 07:29:14  marwan
% help text updated
%
% Revision 3.2  2005/03/16 11:19:02  marwan
% help text modified
%
% Revision 3.1  2004/11/10 07:07:35  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.


error(nargchk(1,1,nargin));
if nargout>2, error('Too many output arguments'), end

warning off
if any(x(:))

  if min(size(x))>100000       % this should speed up the routine; the value
                             % depends on the available memory
    x2=uint8(x);
    N=size(x2);
    x3=zeros(2*N(2)+N(1),N(2));
    x3(N(2)+1:N(2)+N(1),1:N(2))=x2;
    N3=size(x3);
    
    i2=repmat(((1:1+N(2))+N(1)+N(2))',1,N(2));
    i4=i2+repmat((2*N(2)+N(1)+1)*[0:N(2)-1],size(i2,1),1);
    i4(:,end)=[];
    i4=reshape(i4,size(i4,1)*size(i4,2),1);
    x3(i4)=[];
    x3(end)=[];
    x2=(reshape(x3,N(1)+N(2),N(2)))';
  
    x2(end+1,:)=0;
    x=reshape(x2,size(x2,1)*size(x2,2),1);
    x2=x(2:end);x(end)=[];
    z0=find(x==0&x2==1);
    z1=find(x2==0&x==1);
  
  else
  
    N=size(x);
  %   x3=zeros(2*N(2)+N(1),N(2));
  %   x3(N(2)+1:N(2)+N(1),1:N(2))=x;
  %   N3=size(x3);
  %   
  %   i2=repmat(((1:1+N(2))+N(1)+N(2))',1,N(2));
  %   i4=i2+repmat((2*N(2)+N(1)+1)*[0:N(2)-1],size(i2,1),1);
  %   i4(:,end)=[];
  %   i4=reshape(i4,size(i4,1)*size(i4,2),1);
  %   x3(i4)=[];
  %   x3(end)=[];
  %   x=(reshape(x3,N(1)+N(2),N(2)))';
  %  
  %   x(end+1,:)=0;
    
  %  for i1=-ceil(N(2)/2):ceil(N(2)/2); temp=diag(x,i1); X(1:length(temp),1+i1+ceil(N(2)/2))=temp;
  %  end, x=double(X);
    x1=spdiags(x);
    z=reshape(x1,size(x1,1)*size(x1,2),1);
    z2(2:length(z)+1)=z;z2(1)=0;z2(end+1)=0;
    z=diff(z2);
    z0=find(z==1);
    z1=find(z==-1);
  
  end
  if length(z0)>length(z1), z0(end)=[]; end
  if length(z1)>length(z0), z1(end)=[]; end
  
  if isempty(z0), z0=0; end
  if isempty(z1), z1=0; end
  
  if z0(1)>z1(1)
    z0(2:end+1)=z0(1:end);z0(1)=0; 
    if length(z0)>length(z1) 
       z0(end)=[];
    end
  end

  l=sort(z1-z0); %l(end)=[];
  l1=l(find(l-1));
  
  if nargout==2
     b_out=zeros(length(l),1);
     b_out=l';
  end
  
  if nargout>0
     a_out=mean(l1);
  else
     mean(l1)
  end
  
else

  if nargout==2
     b_out=NaN;
  end

  if nargout>0
     a_out=NaN;
  else
     NaN
  end

end

warning on
