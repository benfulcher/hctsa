function varargout=histn(varargin)
%HISTN   N-dimensional histogram.
%    P=HISTN(X) bins the N-dimensional density of X(:,1),...
%    X(:,N) into a 10x10 equally spaced matrix and returns it 
%    in P.
%
%    P=HISTN(X1,X2,...,Xn) bins the M-dimensional density of 
%    X1(i),...,XN(i) into a 10x10 equally spaced matrix and 
%    returns it in P, where M is the sum of the numbers of
%    columns of the input arguments.
%
%    P=HISTN(X1,X2,...,Xn,L), where L is a scalar, uses a lag L
%    between the input vectors. When only one input vector is
%    given, the lag has an effect only if this vector has only 
%    one column. The latter corresponds to P=HISTN(X,X,L).
%
%    P=HISTN(X1,X2,...,Xn,K,L), where K and L are scalars, uses 
%    K bins and a lag L.
%
%    [P,J]=HISTN(...) returns the matrix P and the N-dimensional
%    vector J containing the N-dimensional density matrix and the
%    bin location for X1,...,Xn.
% 
%    HISTN(...) without any output arguments produces a histogram plot.
%
%    Examples: x = randn(10000,3);
%              histn(x)
%
%    See also HIST, HIST2, BAR3, MI.

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 1998-2008
% Andre Sitz, Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:32:57 $
% $Revision: 1.13 $
%
% $Log: histn.m,v $
% Revision 1.13  2009/03/24 08:32:57  marwan
% copyright address changed
%
% Revision 1.12  2007/12/20 16:26:06  marwan
% changed gpl splash behaviour
%
% Revision 1.11  2007/05/25 13:02:24  marwan
% *** empty log message ***
%
% Revision 1.10  2006/10/24 14:16:16  marwan
% minor change: sigma in title line of RP shown only for normalised data
%
% Revision 1.9  2004/11/10 07:05:37  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.



%try 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check the input

error(nargchk(1,16,nargin));
if nargout>2, error('Too many output arguments'), end
nogui=1;
lag=[]; 
nbin=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read the input

x={}; sm_length=[]; sum_m=0;

if isnumeric(varargin{1})
  for i=1:nargin
    if isnumeric(varargin{i}) & max(size(varargin{i}))>1
      y=varargin{i};
      if size(y,1)==1, y=y'; end 
      my=size(y,2); sum_m=sum_m+my;
      mx=size(y,1); 
      if isempty(sm_length), sm_length=mx; end
      if mx < sm_length, sm_length=mx; end
      x(i)={y};
    end
  end
  m=size(x,2); if sum_m==1; lag=1; x(2)=x(1); sum_m=2; m=2; end
  if isempty(x), error('Not a valid input vector.'), end

  i_double=find(cellfun('isclass',varargin,'double'));
  i_char=find(cellfun('isclass',varargin,'char'));
  if ~isempty(i_char)
    if findstr(lower(varargin{i_char(1)}),'n')
      nogui=1; 
    end
    if findstr(lower(varargin{i_char(1)}),'s')
      nogui=2;
    end
    if findstr(lower(varargin{i_char(1)}),'g')
      nogui=0;nogui=1;
    end
  end

if length(i_double)>1
  if max(size(varargin{i_double(end-1)}))==1
    nbin=varargin{i_double(end-1)}; else nbin=10;
  end
end
if max(size(varargin{i_double(end)}))==1
  lag=varargin{i_double(end)}; else if isempty(lag), lag=0; end
end

else
  error('Input must be numeric.')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% estimate histogram

%if size(x,1)<size(x,2), x=x';end
for k=1:m;
    minX = repmat(min(x{k}(1:sm_length,:)),length(x{k}(1:sm_length,:)),1);
    rangeX = repmat(max(x{k}(1:sm_length,:)-repmat(min(x{k}(1:sm_length,:)),length(x{k}(1:sm_length,:)),1)),length(x{k}(1:sm_length,:)),1);
    x1(k) = {(x{k}(1:sm_length,:) - minX) ./ rangeX - eps};
    if m>1 & sm_length/(m-1) <= lag, lag = floor(sm_length / (m-1)) - 1; warning(['Lag is too large. Using ',num2str(lag),' instead.']), end
end

temp=0; kl=0;
for k=1:m
  for l=1:size(x1{k},2)
    temp=temp+fix(x1{k}(1+(lag*(k-1)):(end-(lag*(m-k))),l)*nbin)*nbin^(kl);
    kl=kl+1;
  end
end

p=histc(temp,0:(nbin^sum_m-1))';

p2=eval(['reshape(p',repmat([',nbin'],1,sum_m),')/length(x);']);
 
kl=0;
for k=1:m
  for l=1:size(x1{k},2)
    kl=kl+1;
    minx = min(min(x{k}(:,l))); maxx = max(max(x{k}(:,l)));
    if minx == maxx, minx = minx - floor(nbin/2) - 0.5; maxx = maxx + ceil(nbin/2) - 0.5; end
    bins(kl,:)=minx:(maxx-minx)/(nbin-1):maxx+eps;
  end
end

if nargout==0
  p3=eval(['p2(:,:',repmat([',round(nbin/2)'],1,sum_m-2),')']);
  bar3(bins(2,:),p3)
  dhy=bins(2,2)-bins(2,1);
  set(gca,'xlim',[0 nbin+1-.1],'ylim',[bins(2,1)-dhy bins(2,end)+dhy])
  hx=get(gca,'xtick'); dhx=diff(hx); 
  set(gca,'xtick',[1:dhx(1):nbin+1])
  if dhx>1
    tx=num2str(fix(100*[bins(1,1):(bins(1,end)-bins(1,1))/(length(hx)-1):bins(1,end)]')/100);
  else
    tx=num2str(fix(100*[bins(1,:)]')/100);
  end
  set(gca,'xticklabel',tx)
  xlabel('x'), ylabel('y')
elseif nargout==1
  varargout{1}=p2;
elseif nargout==2
  varargout{1}=p2;
  varargout{2}=bins';
end
