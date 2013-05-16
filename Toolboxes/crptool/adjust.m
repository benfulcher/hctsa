function [x, y]=adjust(a, b, flag)
%ADJUST   Adjusts two two-column vectors.
%    [X, Y]=ADJUST(A, B, FLAG) adjusts the two-column vectors A and B 
%    to the same time scale (in first column)
%
%    FLAG = 0 (default) adjustment by cutting
%           1 adjustment by cubic interpolating
%          -1 adjustment by cubic interpolating and forced length (given by A)
%           2 gap filling by histogram estimation (experimental status)
%           3 gap filling by AR(p) estimation (experimental status)
%
%    Except for FLAG=0, X and Y will have the same length.
%
%    Examples: x = (1:110)';
%              y1 = x(11:end); y1(:,2) = sin(x(11:end) / 10);
%              y2 = x(1:100) / 2; y2(:,2) = sin(x(1:100) / 5);
%              [z1 z2] = adjust(y1,y2);

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 2001-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:30:23 $
% $Revision: 2.3 $
%
% $Log: adjust.m,v $
% Revision 2.3  2009/03/24 08:30:23  marwan
% copyright address changed
%
% Revision 2.2  2005/03/16 11:19:02  marwan
% help text modified
%
% Revision 2.1  2004/11/10 07:07:05  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

error(nargchk(2,3,nargin));
if nargin==2
  flag=0;
end
mflag=0;

Na=size(a); Nb=size(b);

if min(Na)<2 | min(Nb)<2
  error('Time scales have to be specified (first column).')
end

if Na(2)>Na(1)
   a=a';
end

if Nb(2)>Nb(1)
   b=b';
end

if flag~=2 & flag>-1

  [Tmin I2]=max([a(1,1),b(1,1)]);
  if I2==2
    [J1 J2]=min(abs(a(:,1)-Tmin));
    a(1:J2-1,:)=[];
  else
    [J1 J2]=min(abs(b(:,1)-Tmin));
    b(1:J2-1,:)=[];
  end  

  [Tmax I2]=min([a(end,1),b(end,1)]);
  if I2==2
    [J1 J2]=min(abs(a(:,1)-Tmax));
    a(J2+1:end,:)=[];
  else
    [J1 J2]=min(abs(b(:,1)-Tmax));
    b(J2+1:end,:)=[];
  end  
  clear I2 J1 J2

end

switch flag

% cut the unmatched start and end
 case 0
  x=a;
  y=b;

% interpolate missed values to forced length Na
 case -1
  Na=size(a); Nb=size(b);
    x=a;
    if b(end,1)<a(end,1); b(end+1,1)=a(end,1); b(end,2)=b(end-1,2); end
    y(:,1)=x(:,1); 
    for i=0:Nb(2)-2;
       y(:,2+i)=interp1(b(:,1)+.0000001*randn(length(b),1),b(:,2+i),y(:,1),'cubic');
    end

% interpolate missed values
 case 1
  Na=size(a); Nb=size(b);
  if Na(1)>=Nb
    x=a;
    y(:,1)=x(:,1);
    for i=0:Nb(2)-2;
       y(:,2+i)=interp1(b(:,1)+.0000001*randn(length(b),1),b(:,2+i),x(:,1),'cubic');
    end
  else
    y=b;
    x(:,1)=y(:,1);
    for i=0:Na(2)-2;
       x(:,2+i)=interp1(a(:,1)+.0000001*randn(length(a),1),a(:,2+i),y(:,1),'cubic');
    end
  end    

% estimate missed values with histogram
 case 2
   warning off
   a=smoothnan(a);
   b=smoothnan(b);
   [y(:,1),y(:,2)]=filllags_hist(a(:,1),a(:,2),b(:,1),b(:,2));  
   [x(:,1),x(:,2)]=filllags_hist(b(:,1),b(:,2),a(:,1),a(:,2));  
   warning on

% estimate missed values with ar(p)
 case 3
%   warning off
   a=smoothnan(a);
   b=smoothnan(b);
   [y(:,1),y(:,2)]=filllags_ar(a(:,1),a(:,2),b(:,1),b(:,2));  
   [x(:,1),x(:,2)]=filllags_ar(b(:,1),b(:,2),a(:,1),a(:,2));  
   warning on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x0,y0]=filllags_hist(x1,y1,x2,y2)

flag2=2; % 1 - y, 2 - dy

dx=log(diff(x2));
dy=diff(y2);
ni=fix(length(dx)/100);

% berechne die verbundverteilung fuer dx, dy und y
if flag2==1
[pxyz ixyz]=histn(dx,y2(1:end-1),y2(2:end),ni);
else
[pxyz ixyz]=histn(dx,y2(1:end-1),dy,ni);
end

% nur fehlende stuetzstellen anfuegen
A=find(~sum(repmat(x2,1,length(x1))==repmat(x1',length(x2),1)));
% fehlende stuetzstellen anfuegen
x3=[x2; x1(A)];
y3=[y2; repmat(NaN,length(A),1)];

[x3s i]=sort(x3);
y3s=y3(i);

if isnan(y3s(1)), y3s(1)=rand*(max(y3s)-min(y3s))+min(y3s); end
for i=2:length(x3)
  if isnan(y3s(i))
% bestimme dx fuer die fehlstelle
    dx=log(x3s(i)-x3s(i-1));
% bestimme die bin-lokation von dx in der verteilung
    j1=sum(dx>=ixyz(:,1))+1;
    if j1>length(ixyz), j1=length(ixyz); end
% bestimme die bin-lokation von y in der verteilung
    j2=sum(y3s(i-1)>=ixyz(:,2))+1;
    if j2>length(ixyz), j2=length(ixyz); end
    P=permute(pxyz(j1,j2,:),[2,3,1]);
% ersetze die verteilung durch eine gauss-verteilung
    P=P/sum(P);
    a=sum(ixyz(:,3).*P(:));
    b=sqrt(sum((ixyz(:,3)-a).^2.*P(:)));
    ineu=[ixyz(1,3):(ixyz(end,3)-ixyz(1,3))/100:ixyz(end,3)];
    Pneu=normpdf(ineu,a,b);
    [Ps i0]=max(Pneu);
    if flag2==1
       y3s(i)=ineu(i0);
    else
       y3s(i)=y3s(i-1)+ineu(i0);
    end
 end
end

x0=x3s; y0=y3s;
clear A y3s x3s y3 x3 i* j* P* 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x0,y0]=filllags_ar(x1,y1,x2,y2)

ymean=mean(y2);
y2=y2-ymean;
addpath /usr/users4/marwan/matlab/ARfit
A=find(~sum(repmat(x2,1,length(x1))==repmat(x1',length(x2),1)));
% fehlende stuetzstellen anfuegen
x3=[x2; x1(A)];
y3=[y2; repmat(NaN,length(A),1)];
[x3s i]=sort(x3);
y3s=y3(i);
[a b c]=arfit(smoothnan(y2),1,20); m=length(b);
a=[b, c];
disp(num2str(a))
y3s=y3s(:);
a=a(:);

for i=1:m,if isnan(y3s(i)), y3s(i)=rand*(max(y3s)-min(y3s))+min(y3s); end,end
for i=m+1:length(x3)
  if isnan(y3s(i))
    y3s(i)=sum(a(1:end-1).*y3s(i-[1:m]))+a(end)*randn;
 end
end
x0=x3s; y0=y3s+ymean;
