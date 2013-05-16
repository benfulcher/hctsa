function y=trafo(x,a)
%TRAFO   Transforms data to a desired distribution.
%    Y=TRAFO(X,A) transforms the data in vector X to data Y of a desired 
%    distribution Y, where 
%    A=0 normal distribution (default),
%    A=1 uniform distribution, 
%    A=2 exponential distribution
%
%    Example: x = rand(5000,1); 
%             subplot(2,1,1), hist(x,20)   % uniformly distributed
%             y = trafo(x,0); 
%             subplot(2,1,2), hist(y,20)   % normally distributed

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 2001-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:33:47 $
% $Revision: 2.3 $
%
% $Log: trafo.m,v $
% Revision 2.3  2009/03/24 08:33:47  marwan
% copyright address changed
%
% Revision 2.2  2005/03/16 11:19:02  marwan
% help text modified
%
% Revision 2.1  2004/11/10 07:09:23  marwan
% initial import
%
%

error(nargchk(1,2,nargin));
if nargout>1, error('Too many output arguments'), end
if nargin~=2
  a=0;
end

if ~isnumeric(a)
  error('The arguments must be numeric. ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transformation to desired distribution

my=0;
s=1;
x=x(:);
N=size(x);
if N(2)==2
  j=2;
  y(:,1)=x(:,1);
else
  j=1;
end


switch a

  case 0
    i=(-5:10/(N(1)-1):5)';
    fs=cumsum(10*(1/(sqrt(2*pi)*s)).*exp(-(i-my).^2./(2*s^2)));
    verteilung=interp1(fs,i,[1:length(fs)],'nearest')';

  case 1
    verteilung=(-1:2/(N(1)-1):1)';
  
  case 2
    verteilung=sort(exprnd(expfit(hist(x(:,j),30)),N(1),1));

end
[xs i]=sort(x(:,j));
y(i,j)=verteilung;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% replaces NaN with interpolated values 
% (c) amron 2001  -> -> ->  smoothnan.m

N=size(y);

for i=1:N(2),
  y1=y(:,i);
  jnan=find(isnan(y(:,i))); j0=find(~isnan(y(:,i)));
  scale=[1:N(1)]';
  if ~isempty(jnan)
    if jnan(1)==1
       y1=x(j0(1)); y1(2:N(1)+1)=y(:,i); 
       scale=0; scale(2:N(1)+1)=[1:N(1)]';
       jnan=find(isnan(y1));
    end
    if jnan(end)==N(1)
       y1(N(1)+1)=y(j0(end),i);
       scale(N(1)+1)=N(1)+1;
       jnan=find(isnan(y1));
    end
    y1(jnan)=[];
    scale(jnan)=[];
    z=interp1(scale,y1,[1:N(1)]);
    y(:,i)=z';
  end
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalizes the values

y(:,j)=(y(:,j)-mean(y(:,j)))/std(y(:,j));



