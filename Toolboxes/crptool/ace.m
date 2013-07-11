function varargout=ace(x,y,w,ii,oi)
%ACE   Finds optimal transformation and maximal correlation.
%    MCOR=ACE(X,Y [,W,II,OI]) estimates the maximal correlation 
%    of the system theta(X)=phi(Y), where X is one-dimensional
%    and Y can be multi-dimensional.
%
%    [THETA, PHI]=ACE(X,Y [,W,II,OI]) estimates the optimal 
%    transformations THETA and PHI of the system theta(X)=phi(Y).
%
%    [THETA, PHI, MCOR]=ACE(X,Y [,W,II,OI]) estimates the optimal 
%    transformations THETA and PHI and the maximal correlation 
%    MCOR of the system theta(X)=phi(Y).
%
%    [THETA, PHI, MCOR, I, O, Imax, Omax]=ACE(X,Y [,W,II,OI]) 
%    estimates the THETA, PHI and MCOR and outputs the number of
%    inner iterations I, break-up number of inner inner iterations,
%    number of outer iterations O and break-up number of outer
%    inner iterations. If the algorithm doesn't converge, the
%    number of iterations will be negative signed.
%
%    ACE(...) without parameters plots the optimal transformations
%    THETA and PHI.
%
%    Optional parameters:
%    W is the half-length of the boxcar window, II is the maximal 
%    number of inner iterations, OI is the minimal number of outer 
%    iterations. If W=[], the default boxcar window size is 11.
%
%    Examples: x = (-1:.002:1) + .3 * rand(1,1001);
%              y = (-1:.002:1) .^ 2 + .3* rand(1,1001);
%              corrcoef(x,y)
%              ace(y,x)
%
%    References:
%    Breiman, L., Friedman, J. H.:
%    Estimating Optimal Transformations for Multiple regression
%    and Correlation, J. Am. Stat. Assoc., Vol. 80, No. 391, 1985.
%    Voss, H., Kurths, J.: 
%    Reconstruction of nonlinear time delay models from data by the 
%    use of optimal transformations, Phys. Lett. A, 234, 1997.
%
%    See also: MCF

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 2001-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:30:23 $
% $Revision: 1.7 $
%
% $Log: ace.m,v $
% Revision 1.7  2009/03/24 08:30:23  marwan
% copyright address changed
%
% Revision 1.6  2005/03/16 11:19:02  marwan
% help text modified
%
% Revision 1.5  2004/11/10 07:05:02  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check the in/out

error(nargchk(2,5,nargin))
error(nargoutchk(0,7,nargout))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% splash the GPL

filename='ace';
which_res=which([filename,'.m']);
gplrc_path=[strrep(which_res,[filename,'.m'],''), 'private'];
gplrc_file=[gplrc_path, filesep, '.gpl.',filename];
if ~exist(gplrc_path)
  mkdir(strrep(which_res,[filename,'.m'],''),'private')
end
if ~exist(gplrc_file)
  fid=fopen(gplrc_file,'w');
  fprintf(fid,'%s\n','If you delete this file, the GNU Public License will');
  fprintf(fid,'%s','splash up at the next time the programme starts.');
  fclose(fid);

  if exist('gpl')
    txt=gpl;
  else
    txt={'The GNU General Public License was not found.'};
  end
  h=figure('NumberTitle','off',...,
         'ButtonDownFcn','close',...
         'Name','GNU General Public License');
  ha=get(h,'Position');
  h=uicontrol('Style','Listbox',...
            'ButtonDownFcn','close',...
            'CallBack','close',...
            'Position',[0 0 ha(3) ha(4)],...
	    'FontName','Courier',...
	    'BackgroundColor',[.8 .8 .8],...
	    'String',txt);
  waitfor(h)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initialization

flag=1;

if size(y,1)<size(y,2)
  y=y';
end
if size(x,1)<size(x,2)
  x=x';
end

if size(x,2)>1, error('The dimension of x must be one.'); end

if nargin < 5 | isempty(oi)
  oi=1000;
end
if nargin < 4 | isempty(ii)
  ii=100;
end
if nargin < 3 | isempty(w)
  w=round(.1*size(y,1));w=5;
end

N=size(y,1);
dim=size(y,2);

if size(x,1)~=N 
  error('The lengths of x and y must match.') 
end

ocrit=1*eps; icrit=1*eps;      % relative accuracy of the CPU

try 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sorting

[xs ix]=sort(x); [ys iy]=sort(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IID distributed data

xr(ix,1)=(1:N)';
for d=1:dim, yr(iy(:,d),d)=(1:N)'; end

theta=(xr-(N+1)/2)/sqrt(N*(N+1)/12);
phi=(yr-(N+1)/2)/sqrt(N*(N+1)/12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% iteration process

ieps=1.; oeps=1.; 
o=1; ocrit1=1;
maxi=0;

h0=waitbar(0,'','Name','Iteration`s Progress');
pos=get(h0,'Position'); set(h0,'Position',[pos(1) pos(2)+100 pos(3) pos(4)]);
h1=get(get(h0,'Children'),'Title');

while o<=oi & ocrit1>ocrit
%  waitbar(log(eps/ocrit1)/36+1), set(h1,'String',['convergence: ', num2str(ocrit1)])
  waitbar(o/oi), set(h1,'String',['convergence: ', num2str(ocrit1)])
  o=o+1;
  i=1; icrit1=1;

  while i<=ii & icrit1>icrit

    i=i+1; 

    for d=1:dim; sum0=0;
      
      for dd=1:dim 
        if dd~=d sum0=sum0+phi(:,dd); end; 
      end;
      
      A=theta-sum0;
      A=A(iy(:,d));		% A=phi(xk) = sortierung
      ww=-w:w;
%      r=conv(A,exp((-ww.^2)/(.5*(ww(end)-ww(1)))^2));	% faltung = gleitender mittelwert
      r=conv(A,ones(2*w+1,1));	% faltung = gleitender mittelwert
      r=r(w+1:N+w)/(2*w+1);  	% bringt r wieder auf die laenge von y
      
      switch flag
      case 1
         r(1:w)=interp1([w+1:2*w+1],r([w+1:2*w+1]),[1:w],'linear','extrap');
         r(N-w+1:N)=interp1([N-2*w:N-w],r([N-2*w:N-w]),[N-w+1:N],'linear','extrap');
      case 2
         r(1:w)=mean(r(w+1:w+2)); r(N-w+1:N)=mean(r(N-w-1:N-w));
      end

%      phi(:,d)=r(:,d);
      phi(:,d)=r(yr(:,d));	% sortiert r wieder auf ausgangssortierung

    end;
    
    icrit1=ieps;
    if dim==1 sum0=phi; 
    else sum0=sum(phi')'; 
    end;
    ieps=sum((sum0-theta).^2)/N;
    icrit1=abs(icrit1-ieps);	% konvergiert die varianz noch?
      
  end;
  
  A=sum0(ix);
%   A=sum0;
%  r=conv(A,exp((-ww.^2)/(.5*(ww(end)-ww(1)))^2));	% faltung = gleitender mittelwert
  r=conv(A,ones(2*w+1,1));
  r=r(w+1:N+w)/(2*w+1);

  switch flag
  case 1
     r(1:w)=interp1([w+1:2*w+1],r([w+1:2*w+1]),[1:w],'linear','extrap');%
     r(N-w+1:N)=interp1([N-2*w:N-w],r([N-2*w:N-w]),[N-w+1:N],'linear','extrap');
  case 2
     r(1:w)=mean(r(w+1:w+2)); r(N-w+1:N)=mean(r(N-w-1:N-w));
  end

  theta=r(xr);
%   theta=r;
  theta=(theta-mean(theta))/std(theta); 
  ocrit1=oeps; 
  oeps=sum((sum0-theta).^2)/N; 
  ocrit1=abs(ocrit1-oeps);
  
  if maxi<i, maxi=i; end

end  

waitbar(1), if ocrit1>ocrit, set(h1,'String','doesn`t converge!'), ocrit1=-ocrit1; pause(.8), end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output

temp=corrcoef(theta,sum0); 
if nargout==7
  varargout(7)={ocrit1};
end 
if nargout>=6
  varargout(6)={icrit1};
end 
if nargout>=5
  varargout(5)={o-1};
end 
if nargout>=4
  varargout(4)={maxi-1};
end
if nargout>=3
  varargout(3)={temp(1,2)};
end
if nargout>=2;
  varargout(1)={theta};
  varargout(2)={phi};
end
if nargout==1;
  varargout(1)={temp(1,2)};
elseif nargout==0
   a=theta;
   b=phi;
   c=temp(1,2);
   d=maxi-1;
   e=o-1;
   f=icrit1;
   g=ocrit1;
   N=size(b,1);
   dim=size(b,2);

   figure
   subplot(ceil((dim+2)/2),2,1)
   plot(x,a,'k.','MarkerSize',0.5)
   if isempty(inputname(1)); 
      tx=['\Theta(x)']; 
      xlabel('x')
   else 
      tx=['\Theta(',inputname(1) ,')'];
      xlabel(inputname(1))
   end
   ylabel(tx)

   for i=2:dim+1,
      [s1 s2]=sort(y(:,i-1)); 
      thetainv=interp1(a+.0000001.*randn(length(a),1),x,b(s2,i-1),'nearest'); 
      
      subplot(ceil((dim+2)/2),2,i)
      [hAX h1 h2]=plotyy(y(:,i-1),b(:,i-1), s1,thetainv);
      set(h1,'Color','k','LineStyle','none','Marker','.','MarkerSize',0.5)
      set(h2,'Color',[.6 .6 .6],'LineStyle','none','Marker','.','MarkerSize',0.5)
      set(hAX(1),'ycolor','k'), set(hAX(2),'ycolor',.5.*[.6 .6 .6])
      if isempty(inputname(2)); 
         tx=['\Phi_{' num2str(i-1) '}(y)'];
         tx2=['\Theta^{-1}(\Phi_{' num2str(i-1) '}(y))'];
         xlabel('y_*')
      else 
         tx=['\Phi_{' num2str(i-1) '}(',inputname(2) ,')'];
         tx2=['\Theta^{-1}(\Phi_{' num2str(i-1) '}(',inputname(2) ,'))'];
         xlabel(inputname(2))
      end
      ylabel(tx)
      s1(find(isnan(thetainv)))=[];
      thetainv(find(isnan(thetainv)))=[];
      [smax sind]=max(s1);
      h1=text(smax+.051*abs(max(s1)-min(s1)),thetainv(sind),tx2,'FontSize',8);
      set(h1,'Parent',hAX(2))
   end
   sig=0;
   
   h=subplot(ceil((dim+2)/2),2,2*ceil((dim+2)/2));

   txt=char;
   if f<0 | g<0
     txt='doesn''t converge!';
   end
      
   tf=f/eps;tg=g/eps;
   if tf>9999, tf=num2str(tf,'%1.0e'); else, tf=num2str(tf); end
   if tg>9999, tg=num2str(tg,'%1.0e'); else, tg=num2str(tg); end
   
   set(h,'Visible','off')
   text(.1,-.18,date,'FontSize',8,'FontAngle','Italic')
   text(.1,.8,['A C E'],'FontSize',14,'FontWeight','Bold','FontName','new century schoolbook')
   text(.1,.6,['Max. Correlation \Psi: ' num2str(c)])
   text(.1,.5,['5% Significance level: ...' ])   %num2str(sig)])
   text(.1,.35,['Data length: ' num2str(N)])
   text(.1,.25,['Window length: ' num2str(2*w+1)])
   text(.1,.15,['Number of Iterations: ', num2str(d), '/ ',  num2str(e)])
   text(.1,.05,['divergence criteria/eps: ', tf, '/ ',  tg])
   if ~isempty(txt)
      text(.1,-.05,txt,'Color',[.8 0 0],'FontWeight','Bold')
   end
   
end


close(h0)

catch
   delete(h0)
   if ~strcmpi(lasterr,'Interrupt')
     disp('Could not compute the optimal transformations.')
     disp('Try other input arguments.')
   end
   for i=1:nargout,varargout(i)={NaN};end
end
