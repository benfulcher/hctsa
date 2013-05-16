function out=crqad(varargin)
% CRQAD   Computes and plots the diagonalwise CRQA measures.
%    Y=CRQAD(X [,Y] [,param1,param2,...]) 
%    Recurrence quantification analysis of diagonals in the
%    cross recurrence plot of the vectors X and Y as well as
%    X and -Y. The output is a structure (see below).
%
%    The input vectors can be multi-column vectors, where
%    each column will be used as a component of the 
%    phase-space vector. However, if the first column is
%    monotonically increasing, it will be used as an
%    time scale for plotting.
%
%    Y=CRQAD(X,M,T,E,W) computes the recurrence
%    quantification analysis of the recurrence plot
%    of X by using the dimension M, delay T, the
%    size of neighbourhood E, for the diagonals within
%    the range [-W,W] around the main diagonal.
%
%    Parameters: dimension M, delay T, the size of
%    neighbourhood E and the range size W are the first 
%    five numbers after the data series; if W is empty,
%    the whole plot will be calculated. Further parameters 
%    can be used to switch between various methods of finding
%    the neighbours of the phasespace trajectory, to suppress 
%    the normalization of the data and to suppress the GUI
%    (useful in order to use this programme by other programmes).
%
%    Methods of finding the neighbours.
%      maxnorm     - Maximum norm.
%      euclidean   - Euclidean norm.
%      minnorm     - Minimum norm.
%      nrmnorm     - Euclidean norm between normalized vectors
%                    (all vectors have the length one).
%      rr          - Maximum norm, fixed recurrence rate.
%      fan         - Fixed amount of nearest neighbours.
%      inter       - Interdependent neighbours.
%      omatrix     - Order matrix.
%      opattern    - Order patterns recurrence plot.
%
%    Normalization of the data series.
%      normalize   - Normalization of the data.
%      nonormalize - No normalization of the data.
%
%    Parameters not needed to be specified.
%
%    Output:
%      Y.RRp  = RRp
%      Y.RRm  = RRm
%      Y.DETp = DETp
%      Y.DETm = DETm
%      Y.Lp   = Lp
%      Y.Lm   = Lm
%      
%    Examples: a = sin(0:.1:80) + randn(1,801);
%              b = sin(0:.1:80) + randn(1,801);
%              crqad(a,b,3,15,.1,100,'fan')
%
%    See also CRQA, CRQAD_BIG, CRP, CRP2, CRP_BIG, DL, TT, RPDE.
%
%    References: 
%    Marwan, N., Kurths, J.:
%    Nonlinear analysis of bivariate data with cross recurrence plots,
%    Phys. Lett. A, 302, 2002.

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 2002-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2010/06/30 12:03:02 $
% $Revision: 2.10 $
%
% $Log: crqad.m,v $
% Revision 2.10  2010/06/30 12:03:02  marwan
% Help text modified
%
% Revision 2.9  2009/03/24 08:35:19  marwan
% corrected XCF calculation
%
% Revision 2.8  2007/07/18 17:18:44  marwan
% integer values in the arguments supported
%
% Revision 2.7  2007/05/15 17:33:13  marwan
% new neighbourhood criterion: fixed RR
%
% Revision 2.6  2006/10/24 14:16:16  marwan
% minor change: sigma in title line of RP shown only for normalised data
%
% Revision 2.5  2006/07/04 14:03:57  marwan
% axis-error
%
% Revision 2.4  2005/03/16 11:19:02  marwan
% help text modified
%
% Revision 2.3  2004/11/12 08:40:46  marwan
% order patterns recurrence plot added
%
% Revision 2.2  2004/11/10 07:05:55  marwan
% initial import
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% programme properties

global props

init_properties

lmin=10;
w=[]; method='max'; method_n=1; t=1; m=1; e=.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check the input

error(nargchk(1,10,nargin));
if nargout>1, error('Too many output arguments'), end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% splash the GPL

splash_gpl('crp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check and read the input

   varargin{11}=[];
   % transform any int to double
   intclasses = {'uint8';'uint16';'uint32';'uint64';'int8';'int16';'int32';'int64'};
   flagClass = [];
   for i = 1:length(intclasses)
       i_int=find(cellfun('isclass',varargin,intclasses{i}));
       if ~isempty(i_int)
           for j = 1:length(i_int)
               varargin{i_int(j)} = double(varargin{i_int(j)});
           end
           flagClass = [flagClass; i_int(:)];
       end
   end
   if ~isempty(flagClass)
       disp(['Warning: Input arguments at position [',num2str(flagClass'),'] contain integer values']);
       disp(['(now converted to double).'])
   end
   i_double=find(cellfun('isclass',varargin,'double'));
   i_char=find(cellfun('isclass',varargin,'char'));
   nogui=0;

   if nargin & isnumeric(varargin{1})

     % check the text input parameters for method, gui 
      check_meth={'ma','eu','mi','nr','rr','fa','in','om','op','di'}; 	% maxnorm, euclidean, nrmnorm,  fan, distance
      check_gui={'gui','nog','sil'};		% gui, nogui, silent
      temp_meth=0;
      temp_gui=0;
      if ~isempty(i_char)
         for i=1:length(i_char), 
            varargin{i_char(i)}(4)='0';
            temp_gui=temp_gui+strcmpi(varargin{i_char(i)}(1:3),check_gui'); 
            temp_meth=temp_meth+strcmpi(varargin{i_char(i)}(1:2),check_meth'); 
         end
         method_n=min(find(temp_meth));
         nogui=min(find(temp_gui))-1;
         for i=1:length(i_char); temp2(i,:)=varargin{i_char(i)}(1:3); end
         i_char(strmatch(check_gui(find(temp_gui)),temp2))=[];
         if isempty(nogui), nogui=0; end
         if isempty(method_n), method_n=1; end
         if nogui>2, nogui=1; end
         if method_n>length(check_meth), method0=length(check_meth); end
         method=check_meth{method_n};
      else
         nogui=0;
         if nargout
            nogui=1;
            action='compute'; 
         end
      end
      if nogui==0
        action='init';
      else
        action='compute'; 
      end
      
      % get the parameters for creating RP
      if max(size(varargin{1}))<=3
         error('To less values in data X.')
      end
      x=double(varargin{1});
      if isempty(varargin{2}) | ~isnumeric(varargin{2}), y=x; else
      y=double(varargin{2}); end
      if sum(double(diff(x(:,1))<=0)), embed_flag=0; end
   
      if (isnumeric(varargin{2}) & max(size(varargin{2}))==1) | ~isnumeric(varargin{2})
        y=x;
        if ~isempty(varargin{i_double(2)}), m=varargin{i_double(2)}(1); else m=1; end
        if ~isempty(varargin{i_double(3)}), t=varargin{i_double(3)}(1); else t=1; end
        if ~isempty(varargin{i_double(4)}), e=varargin{i_double(4)}(1); else e=.1; end
        if ~isempty(varargin{i_double(5)}), w=varargin{i_double(5)}(1); else w=varargin{i_double(5)}; end
%        if ~isempty(varargin{i_double(6)}), wstep=varargin{i_double(6)}(1); else wstep=1; end
      else
        if ~isempty(varargin{i_double(3)}), m=varargin{i_double(3)}(1); else m=1; end
        if ~isempty(varargin{i_double(4)}), t=varargin{i_double(4)}(1); else t=1; end
        if ~isempty(varargin{i_double(5)}), e=varargin{i_double(5)}(1); else e=.1; end
        if ~isempty(varargin{i_double(6)}), w=varargin{i_double(6)}(1); else w=varargin{i_double(6)}; end
%        if ~isempty(varargin{i_double(7)}), wstep=varargin{i_double(7)}(1); else wstep=1; end
      end
    else
      error('No valid arguments.')
    end

   
    Nx=length(x); Ny=length(y);
    if size(x,1)<size(x,2), x=x'; end
    if size(y,1)<size(y,2), y=y'; end
   
   if size(x,2)>=2
      xscale=x(:,1); 
      if ~isempty(find(diff(xscale)<0)), embed_flag=0;end
   else
      xscale=(1:length(x))'; 
   end
   if size(y,2)>=2
      yscale=y(:,1); 
      if ~isempty(find(diff(yscale)<0)), embed_flag=0;end
   else
       yscale=(1:length(y))';
   end
      
      if max(size(x))~=max(size(y)),
        if ~nogui, errordlg('Data must have the same length.','Check Data'), else error('Data must have the same length.'), end
      end
      if e<0, 
        e=1; 
	if ~nogui
	   warndlg('The threshold size E can not be negative and is now set to 1.','Check Data')
	   h=findobj('Tag','crqa_eps');
	   set(h(1),'String',str2num(e)) 
        else 
	   disp('The threshold size E can not be negative and is now set to 1.'), 
        end
      end
      if t<1, 
        t=1; 
	if ~nogui
	   warndlg('The delay T can not be smaller than one and is now set to 1.','Check Data')
	   h=findobj('Tag','crqa_maxLag');
	   set(h(1),'String',str2num(t)) 
	else
	   disp('The delay T can not be smaller than one and is now set to 1.')
	end
      end
      if isempty(w), w=.5*Nx; wstep=1; end
%      if w<2, 
%        w=2; 
%	if ~nogui, warndlg('The window size W exceeds the valid range.','Check Data')
%	else, disp('The window size W exceeds the valid range.'), end
%      end
      if w>Nx, 
        w=Nx; wstep=1;; 
	if ~nogui, warndlg('The window size W exceeds the valid range.','Check Data')
	else, disp('The window size W exceeds the valid range.'), end
      end
      t=round(t); m=round(m); w=round(w);% wstep=round(wstep); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute

flag=1;

x1=x; x2=y;
  if length(method)>1 & strcmpi(method(1:2),'di')
      disp('Warning: RQA from distance plot not possible!')
      return
  end

warning off
if size(x1,1)<size(x1,2), x1=x1'; end
if size(x2,1)<size(x2,2), x2=x2'; end

x=crp2(x1,x2,m,t,e,method,'sil','nor');
warning off
  N=size(x);
  x3=zeros(2*N(2)+N(1),N(2));
  x3(N(2)+1:N(2)+N(1),1:N(2))=x;
  N3=size(x3);
  
  i2=repmat(((1:1+N(2))+N(1)+N(2))',1,N(2));
  i4=i2+repmat((2*N(2)+N(1)+1)*[0:N(2)-1],size(i2,1),1);
  i4(:,end)=[]; i4=reshape(i4,size(i4,1)*size(i4,2),1);
  x3(i4)=[]; x3(end)=[];
  x4=(reshape(x3,N(1)+N(2),N(2)))'; x4(end+1,:)=0;

i=1; clear DET L RR
for j=-w:w, clear z z0 z1

  if sum(x4(:,N(2)+1+j))==N(2), l1=N(2);else
  z=diff(x4(:,N(2)+1+j)); 
  if ~isempty(find(~(z-1))),z0(:,1)=find(~(z-1));else,z0(1:N(2))=0;end, 
  if ~isempty(find(~(z+1))),z1=find(~(z+1));else,z1(1:N(2))=0;end
  
  if z0(1)>z1(1)
    z0(2:end+1,1)=z0(1:end);z0(1)=0; 
    if length(z0)>length(z1), z0(end)=[]; end
  end
  l=sort(z1-z0); l1=l(find(l>lmin));
  end
  DET(i)=sum(l1)/sum(x4(:,N(2)+1+j));
  L(i)=mean(l1);
  RR(i)=sum(x4(:,N(2)+1+j))/(N(2)-abs(j));
  i=i+1;
end
L(find(isnan(L)))=0;
RR(find(isnan(RR)))=0;
DET(find(isnan(DET)))=0;

if nargout, XCF=xcf(x1,x2,w,1); end

if ~nargout
subplot(2,2,1)
clim=1;
xcf(x1(:,1),x2(:,1),w)
set(gca,'fonta','i')
xlabel('Lag'), axis([-w w -clim clim])
ylabel('Cross Correlation')
h=text(0,0,'A','fontw','b');set(h,'un','pi'),set(h,'pos',[9,16,0])

switch flag
case 1
subplot(2,2,2), plot([-w:w],RR,'k','linew',.7), 
set(gca,'fonta','i'),axis([-w w 0 1])
xlabel('Lag'),ylabel('Recurrence Rate'),grid on
h=text(0,0,'B','fontw','b');set(h,'un','pi'),set(h,'pos',[9,16,0])

subplot(2,2,3), plot([-w:w],DET,'k','linew',.7), 
set(gca,'fonta','i'),xlabel('Lag')
axis([-w w 0 1]),ylabel('Determinism'),grid on
h=text(0,0,'C','fontw','b');set(h,'un','pi'),set(h,'pos',[9,16,0])

subplot(2,2,4), plot([-w:w],L,'k','linew',.7), 
set(gca,'fonta','i'),xlabel('Lag')

axis([-w w 0 max([max(L) 1])]),ylabel('Averaged Line Length'),grid on
h=text(0,0,'D','fontw','b');set(h,'un','pi'),set(h,'pos',[9,16,0])

case 2
subplot(2,2,2), plot([-w:w],smooth(RR,5,5),'k','linew',.7), 
set(gca,'fonta','i'),axis([-w w 0 1])
xlabel('Lag'),ylabel('Recurrence Rate'),grid on
h=text(0,0,'B','fontw','b');set(h,'un','pi'),set(h,'pos',[9,16,0])

subplot(2,2,3), plot([-w:w],smooth(DET,5,5),'k','linew',.7), 
set(gca,'fonta','i'),xlabel('Lag')
axis([-w w 0 1]),ylabel('Determinism'),grid on
h=text(0,0,'C','fontw','b');set(h,'un','pi'),set(h,'pos',[9,16,0])

subplot(2,2,4), plot([-w:w],smooth(L,5,5),'k','linew',.7), 
set(gca,'fonta','i'),xlabel('Lag')
axis([-w w 0 max(L)]),ylabel('Averaged Line Length'),grid on
h=text(0,0,'D','fontw','b');set(h,'un','pi'),set(h,'pos',[9,16,0])

end
else
  out.XCF=XCF';
  out.RRp=RR;
  out.DETp=DET;
  out.Lp=L;
end

x=crp2(x1,-x2,m,t,e,method,'sil');
warning off

  N=size(x);
  x3=zeros(2*N(2)+N(1),N(2));
  x3(N(2)+1:N(2)+N(1),1:N(2))=x;
  N3=size(x3);
  
  i2=repmat(((1:1+N(2))+N(1)+N(2))',1,N(2));
  i4=i2+repmat((2*N(2)+N(1)+1)*[0:N(2)-1],size(i2,1),1);
  i4(:,end)=[]; i4=reshape(i4,size(i4,1)*size(i4,2),1);
  x3(i4)=[]; x3(end)=[];
  x4=(reshape(x3,N(1)+N(2),N(2)))'; x4(end+1,:)=0;

i=1; clear DET L RR
for j=-w:w, clear z z0 z1

  if sum(x4(:,N(2)+1+j))==N(2), l1=N(2);else
  z=diff(x4(:,N(2)+1+j)); 
  if ~isempty(find(~(z-1))),z0(:,1)=find(~(z-1));else,z0(1:N(2))=0;end, 
  if ~isempty(find(~(z+1))),z1=find(~(z+1));else,z1(1:N(2))=0;end

  if z0(1)>z1(1)
    z0(2:end+1,1)=z0(1:end);z0(1)=0; 
    if length(z0)>length(z1), z0(end)=[]; end
  end
  l=sort(z1-z0); l1=l(find(l>lmin));
  end
  DET(i)=sum(l1)/sum(x4(:,N(2)+1+j));
  L(i)=mean(l1);
  RR(i)=sum(x4(:,N(2)+1+j))/(N(2)-abs(j));
  i=i+1;
end
L(find(isnan(L)))=0;
RR(find(isnan(RR)))=0;
DET(find(isnan(DET)))=0;

if ~nargout
switch flag
case 1
subplot(2,2,2), hold on, plot([-w:w],RR,'r','linew',.7), hold off, set(gca,'YLimMode','Auto')
subplot(2,2,3), hold on, plot([-w:w],DET,'r','linew',.7), hold off, hold off, set(gca,'YLimMode','Auto')
subplot(2,2,4), hold on, plot([-w:w],L,'r','linew',.7), hold off, hold off, set(gca,'YLimMode','Auto')
case 2
subplot(2,2,2), hold on, plot([-w:w],smooth(RR,5,5),'r','linew',.7), hold off, hold off, set(gca,'YLimMode','Auto')
subplot(2,2,3), hold on, plot([-w:w],smooth(DET,5,5),'r','linew',.7), hold off, hold off, set(gca,'YLimMode','Auto')
subplot(2,2,4), hold on, plot([-w:w],smooth(L,5,5),'r','linew',.7), hold off, hold off, set(gca,'YLimMode','Auto')
end
else
  out.RRm=RR;
  out.DETm=DET;
  out.Lm=L;
end
warning on

