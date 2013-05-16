function xout=erqa(varargin)
%ERQA   Computes and plots the ERQA measures.
%    Y=ERQA(X [,Y] [,param1,param2,...]) 
%    Extended recurrence quantification analysis of 
%    cross-recurrence plots with the first vector 
%    X and the second Y. 
%
%    Y=ERQA(X,M,T,E,W,WS) computes the extended
%    recurrence quantification analysis of the recurrence 
%    plot of X by using the dimension M, delay T, the
%    size of neighbourhood E, the window size W and 
%    a window shifting value of WS.
% 
%    Parameters: dimension M, delay T, the size of
%    neighbourhood E, the window size W and the shift
%    value WS are the first five numbers after the data 
%    series; if W is empty, the whole plot will be calculated.
%    Further parameters can be used to switch between various 
%    methods of finding the neighbours of the phasespace 
%    trajectory, to suppress the normalization of the data 
%    and to suppress the GUI (useful in order to use this 
%    programme by other programmes).
%
%    Methods of finding the neighbours.
%      maxnorm     - Maximum norm.
%      euclidean   - Euclidean norm.
%      nrmnorm     - Euclidean norm between normalized vectors
%                    (all vectors have the length one).
%      fan         - Fixed amount of nearest neighbours.
%      inter       - Interdependent neighbours.
%      distance    - Distance coded matrix (global CRP).
%
%    Normalization of the data series.
%      normalize   - Normalization of the data.
%      nonormalize - No normalization of the data.
%
%    Suppressing the GUI.
%      gui         - Creates (the GUI and) the output plot.
%      nogui       - Suppresses the GUI and the output plot.
%      silent      - Suppresses all output.
%
%    Parameters not needed to specify.
%
%    Output:
%      Y(:,1)=K1
%      Y(:,2)=K2
%      Y(:,3)=D2
%      Y(:,4)=I2
%      
%    Examples: a=randn(300,1);
%              erqa(a,1,1,.2,40,2,'max')
%
%              N=300; w=40; ws=2;
%              a=3.4:.6/(N-1):4;
%              b=.5; for i=2:N, b(i)=a(i)*b(i-1)*(1-b(i-1));end
%              y=erqa(b,3,2,.1,w,ws,'nogui');
%              subplot(2,1,1), plot(a,b,'.','markersize',.1)
%              title('logistic map'), axis([3.4 4 0 1])
%              subplot(2,1,2), plot(a(1:ws:N-w),y(1:ws:N-w,1))
%              ylabel('recurrence rate'), axis([3.4 4 0 1])
%
%      
%    See also CRP, CRP2, CRP_BIG, DL, TT.
%
%    References: 
%    Thiel, M., Romano, M. C., Kurths, J.: 
%    Analytical description of Recurrence Plots of stochastic
%    and chaotic processes, to be publ. 2002

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 1998-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2004/11/10 07:06:03 $
% $Revision: 1.9 $
%
% $Log: erqa.m,v $
% Revision 1.9  2004/11/10 07:06:03  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check and read the input

error(nargchk(1,10,nargin));
if nargout>1, error('Too many output arguments'), end

   varargin{11}=[];
   i_double=find(cellfun('isclass',varargin,'double'));
   i_char=find(cellfun('isclass',varargin,'char'));

   % check the text input parameters for gui 
   check_gui={'gui','nog','sil'};		% gui, nogui, silent
   temp_gui=0;
   if ~isempty(i_char)
      for i=1:length(i_char), 
         varargin{i_char(i)}(4)='0';
         temp_gui=temp_gui+strcmpi(varargin{i_char(i)}(1:3),check_gui'); 
      end
      nogui=min(find(temp_gui))-1;
      for i=1:length(i_char); temp2(i,:)=varargin{i_char(i)}(1:3); end
      i_char(strmatch(check_gui(find(temp_gui)),temp2))=[];
      if isempty(nogui), nogui=1; end
      if nogui>2, nogui=1; end
   else
      nogui=1;
   end
   if nogui==0, warning('Sorry. GUI not yet available.'),end

   % get the parameters for creating RP
   if max(size(varargin{1}))<=3
      error('To less values in data X.')
   end
   x=double(varargin{1});
   if isempty(varargin{2}) | ~isnumeric(varargin{2}), y=x; else
   y=double(varargin{2}); end

   if (isnumeric(varargin{2}) & max(size(varargin{2}))==1) | ~isnumeric(varargin{2})
     y=x;
     if ~isempty(varargin{i_double(2)}), m=varargin{i_double(2)}(1); else m=1; end
     if ~isempty(varargin{i_double(3)}), t=varargin{i_double(3)}(1); else t=1; end
     if ~isempty(varargin{i_double(4)}), e=varargin{i_double(4)}(1); else e=.1; end
     if ~isempty(varargin{i_double(5)}), w=varargin{i_double(5)}(1); else w=varargin{i_double(5)}; end
     if ~isempty(varargin{i_double(6)}), wstep=varargin{i_double(6)}(1); else wstep=1; end
   else
     if ~isempty(varargin{i_double(3)}), m=varargin{i_double(3)}(1); else m=1; end
     if ~isempty(varargin{i_double(4)}), t=varargin{i_double(4)}(1); else t=1; end
     if ~isempty(varargin{i_double(5)}), e=varargin{i_double(5)}(1); else e=.1; end
     if ~isempty(varargin{i_double(6)}), w=varargin{i_double(6)}(1); else w=varargin{i_double(6)}; end
     if ~isempty(varargin{i_double(7)}), wstep=varargin{i_double(7)}(1); else wstep=1; end
   end

   Nx=length(x); Ny=length(y);
   if size(x,1)<size(x,2), x=x'; end
   if size(y,1)<size(y,2), y=y'; end

  if size(x,2)>=2
     xscale=x(:,1); 
     if ~isempty(find(diff(xscale)<0)), error('First column of the first vector must be monotonically non-decreasing.'),end
     x=x(:,2);
  else
     xscale=(1:length(x))'; 
  end
  if size(y,2)>=2
     yscale=y(:,1); 
     if ~isempty(find(diff(yscale)<0)), error('First column of the second vector must be monotonically non-decreasing.'),end
     y=y(:,2);
  else
      yscale=(1:length(y))';
  end
   
   if max(size(x))~=max(size(y)),error('Data must have the same length.'),end
   if e<0, e=1; disp('The threshold size E can not be negative and is now set to 1.'), end
   if t<1, t=1; disp('The delay T can not be smaller than one and is now set to 1.'), end
   if isempty(w), w=Nx-1; end
   if w<2, w=2; disp('The window size W exceeds the valid range.'), end
   if w>Nx, w=Nx-1; disp('The window size W exceeds the valid range.'), end
   if wstep<2 & wstep>Nx/3, wstep=2; disp('The window shifting value WS exceeds the valid range.'), end
   t=round(t); m=round(m); w=round(w); wstep=round(wstep); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% splash the GPL

filename='crp';
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute

%try 

de=.1;
if nogui~=2, h=waitbar(0,'1');h1=get(h,'chil');h1=get(h1,'title'); end
for i=1:wstep:Nx-w; 
if nogui~=2, waitbar(i/(Nx-w)); 
txt=['CRP1 - ',num2str(i),'/',num2str(Nx-w)];
set(h1,'str',txt); drawnow, end

try
  X1=crp(x(i:i+w-1,:),x(i:i+w-1,:),m,t,e,varargin{i_char},'silent');
  if nogui~=2, txt=['CRP1a - ',num2str(i),'/',num2str(Nx-w)]; set(h1,'str',txt); drawnow, end
  X1a=crp(y(i:i+w-1,:),y(i:i+w-1,:),m,t,e,varargin{i_char},'silent');
  X2 = X1 .* X1a; X1 = X2;
  if nogui~=2, txt=['CRP2 - ',num2str(i),'/',num2str(Nx-w)]; set(h1,'str',txt); drawnow, end
 X2=crp(x(i:i+w-1,:),x(i:i+w-1,:),m,t,e+de,varargin{i_char},'silent');
 X1a=crp(y(i:i+w-1,:),y(i:i+w-1,:),m,t,e+de,varargin{i_char},'silent');
  X2a = X2 .* X1a; X2 = X2a;
  if nogui~=2, txt=['CRP3 - ',num2str(i),'/',num2str(Nx-w)]; set(h1,'str',txt); drawnow, end
  X3=crp(x(i:i+w-1,:),x(i:i+w-1,:),m+2,t,e,varargin{i_char},'silent');
  X1a=crp(y(i:i+w-1,:),y(i:i+w-1,:),m+2,t,e,varargin{i_char},'silent');
  X2a = X3 .* X1a; X3 = X2a;
catch
  if nogui~=2, close(h), end 
  error(lasterr)
end

txt=['ERQA - ',num2str(i),'/',num2str(Nx-w)];
if nogui~=2, set(h1,'str',txt); drawnow, end
warning off 

N=length(X1);
[ml NL]=dl(X1);
pl=hist(NL,[0:N]);
for j=1:N;
  plc(j)=sum(pl([1:N-j]+j).*([1:N-j]+1));
end
Np=min(find(plc==0));

% 1st measure
i1=[round(.04*Np):round(.12*Np)]+1; 
i2=[round(.3*Np):round(.85*Np)]+1; 
%i1=2:round(.12*Np); 
if length(i1)==1, i1=[]; end
if length(i2)==1, i2=[]; end
if ~isempty(i1),K1=-regress(log(plc(i1))',[i1',ones(length(i1),1)]);else K1=NaN; end
if ~isempty(i2),K2=-regress(log(plc(i2))',[i2',ones(length(i2),1)]);else K2=NaN; end
%K2=regress(log(plc(i2))',[i2',ones(length(i1),1)]);
%K2=log(plc(i2));

% 2nd measure
txt=['2RQA - ',num2str(i),'/',num2str(Nx-w)];
if nogui~=2, set(h1,'str',txt); drawnow, end
[ml NL]=dl(X2); N2=length(X2);
pl=hist(NL,[0:N2]);
for j=1:N2; plc2(j)=sum(pl([1:N2-j]+j).*([1:N2-j]+1)); end
l=round(ml);
D2=log(plc(l)/plc2(l))/log(e/(e+de));

% 3rd measure
H2=-log(sum(X1(:))/(N^2));
X1=double(X1);
for t=1:3;
  j=1:N-t;k=1:N-t;
  H2t(t)=-log(sum(sum(X1(k,j).*X1(k+t,j+t)))/(N^2));
end

I2=2*H2-H2t;

% ApEn
B=(sum(X1,1)-1)/(length(X1)+1); A=(sum(X3,1)-1)/(length(X3)+1);

ApEn = sum(log(B(find(B))))/(length(X1)+1) - sum(log(A(find(A))))/(length(X3)+1);

%SampEn
B=(sum(X1,1)-1)/(length(X1)-1); A=(sum(X3,1)-1)/(length(X1)-1);
SampEn = -log(mean(A)/mean(B));


warning on

Y(i,1)=K1(1); 
Y(i,2)=K2(1);
Y(i,3)=D2;
Y(i,4)=I2(1);
Y(i,5)=I2(2);
Y(i,6)=I2(3);
Y(i,7)=ApEn;
Y(i,8)=SampEn;

end

if nogui~=2
  close(h)
  tx={'K1';'K2';'D2';'I2';'I2';'I2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot

   for i=1:4,subplot(2,2,i),plot(Y(1:wstep:end,i)),ylabel(tx(i));end
   hold on
   plot(Y(1:wstep:end,i+1),'r')
   plot(Y(1:wstep:end,i+2),'g')
end

if nargout, xout=Y; end

if 0
%catch
  if ishandle(h),close(h),end
  error(lasterr)
end
