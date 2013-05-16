function out=crqad_big(varargin)
% CRQAD_BIG   Computes and plots the diagonalwise CRQA measures.
%    Y=CRQAD_BIG(X [,Y] [,param1,param2,...]) 
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
%    Examples: a = sin(0:.1:800) + randn(1,8001);
%              b = sin(0:.1:800) + randn(1,8001);
%              crqad_big(a,b,3,15,.1,100,'euc')
%
%    See also CRQA, CRQAD, CRP, CRP2, CRP_BIG, DL, TT.
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
%
% $Date: 2009/03/24 08:32:09 $
% $Revision: 5.4 $
%
% $Log: crqad_big.m,v $
% Revision 5.4  2009/03/24 08:32:09  marwan
% copyright address changed
%
% Revision 5.3  2007/07/18 17:18:44  marwan
% integer values in the arguments supported
%
% Revision 5.2  2006/10/24 14:36:06  marwan
% minor bug: different lengths of vectors z0 and z1
%
% Revision 5.1  2006/10/09 07:49:03  marwan
% initial check in
%
%
%
% This program is part of the new generation XXII series.
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% programme properties

global props

init_properties

lmin=3;
w=[]; method='max'; method=1; t=1; m=1; e=.1;

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
      check_meth={'ma','eu','mi','nr','fa','in','om','op','di'}; 	% maxnorm, euclidean, nrmnorm,  fan, distance
      check_gui={'gui','nog','sil'};		% gui, nogui, silent
      check_norm={'non','nor'};				% nonormalize, normalize
      temp_meth=0;
      temp_gui=0;
      temp_norm=0;
      if ~isempty(i_char)
         for i=1:length(i_char), 
            varargin{i_char(i)}(4)='0';
            temp_meth=temp_meth+strcmpi(varargin{i_char(i)}(1:2),check_meth'); 
            temp_gui=temp_gui+strcmpi(varargin{i_char(i)}(1:3),check_gui'); 
            temp_norm=temp_norm+strcmpi(varargin{i_char(i)}(1:3),check_norm'); 
         end
         method=min(find(temp_meth));
         nonorm=min(find(temp_norm))-1;
         nogui=min(find(temp_gui))-1;
         for i=1:length(i_char); temp2(i,:)=varargin{i_char(i)}(1:3); end
         i_char(strmatch(check_gui(find(temp_gui)),temp2))=[];
         if isempty(nogui), nogui=0; end
         if isempty(method), method=1; end
         if nonorm>1, nonorm=1; end
         if nogui>2, nogui=1; end
         if method>length(check_meth), method=length(check_meth); end
      else
         nogui=0; nonorm=1; 
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
      if nonorm==1, x=(x(:,2)-mean(x(:,2)))/std(x(:,2)); else x=x(:,2); end
   else
      if nonorm==1, x=(x-mean(x))/std(x); end
      xscale=(1:length(x))'; 
   end
   if size(y,2)>=2
      yscale=y(:,1); 
      if ~isempty(find(diff(yscale)<0)), embed_flag=0;end
      if nonorm==1, y=(y(:,2)-mean(y(:,2)))/std(y(:,2)); else y=y(:,2); end
   else
       if nonorm==1, y=(y-mean(y))/std(y); end
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
      if isempty(w) & Nx > 5000, w = 100; wstep=1; end
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
  if method > 3
      disp('Warning: RQA from choosen neighbourhood is not possible!')
      return
  end

warning off
if size(x1,1)<size(x1,2), x1=x1'; end
if size(x2,1)<size(x2,2), x2=x2'; end


  
  % embedding vectors
  NX=Nx-t*(m-1);NY=Ny-t*(m-1);  
  i=(1:NX)';j=0:t:0+(m-1)*t;
  i=reshape(i(:,ones(1,m))+j(ones(NX,1),:),m*NX,1);
  x1=x(i);
  x2=reshape(x1,NX,m);

  i=(1:NY)';j=0:t:0+(m-1)*t;
  i=reshape(i(:,ones(1,m))+j(ones(NY,1),:),m*NY,1);
  y1=y(i);
  y2=reshape(y1,NY,m);


% compute the diagonalwise RQA
clear DET L RR; j = 0;
hw = waitbar(0,'Compute RQA');
for i = 0:w,waitbar(i/(2*w))
clear z z0 z1

    if m > 1

        switch(method)

        %%%%%%%%%%%%%%%%% local CRP, fixed distance

          case 1
        %%%%%%%%%%%%%%%%% maximum norm
            s = max(abs(x2(1:NX-i,:) - y2(i+1:NY,:))');
          case 2
        %%%%%%%%%%%%%%%%% euclidean norm
            errcode=112;
            s = sqrt(sum((x2(1:NX-i,:) - y2(i+1:NY,:)).^2, 2));
          case 3
        %%%%%%%%%%%%%%%%% minimum norm
            errcode=113;
            s = sum(abs(x2(1:NX-i,:) - y2(i+1:NY,:))');
        end
    else
        s = abs(x2(1:NX-i,:) - y2(i+1:NY,:));
    end
    
    X = s(:) < e;
    

    if sum(X) == length(X), l1=length(X); else
        z=diff(X); z0 = [];
        if ~isempty(find(~(z-1))),z0(:,1)=find(~(z-1));else,z0(1:length(X))=0;end, 
        if ~isempty(find(~(z+1))),z1=find(~(z+1));else,z1(1:length(X))=0;end

        if z0(1)>z1(1)
            z0(2:end+1,1)=z0(1:end);z0(1)=0; 
        end
        if length(z0)>length(z1), z0(end)=[]; end
        l=sort(z1-z0); l1=l(find(l>lmin));
    end
    DET(j+w+1)=sum(l1)/sum(X);
    L(j+w+1)=mean(l1);
    RR(j+w+1)=sum(X)/length(X);

    if m > 1

        switch(method)

        %%%%%%%%%%%%%%%%% local CRP, fixed distance

          case 1
        %%%%%%%%%%%%%%%%% maximum norm
            s = max(abs(x2(i+1:NX,:) - y2(1:NY-i,:))');
          case 2
        %%%%%%%%%%%%%%%%% euclidean norm
            errcode=112;
            s = sqrt(sum((x2(i+1:NX,:) - y2(1:NY-i,:)).^2, 2));
          case 3
        %%%%%%%%%%%%%%%%% minimum norm
            errcode=113;
            s = sum(abs(x2(i+1:NX,:) - y2(1:NY-i,:))');
        end
    else
        s = abs(x2(i+1:NX,:) - y2(1:NY-i,:));
    end
    X = s(:) < e;

    if sum(X) == length(X), l1=length(X);else
        z=diff(X); z0 = [];
        if ~isempty(find(~(z-1))),z0(:,1)=find(~(z-1));else,z0(1:length(X))=0;end, 
        if ~isempty(find(~(z+1))),z1=find(~(z+1));else,z1(1:length(X))=0;end
            z1 = z1(:); z0 = z0(:);

        if z0(1)>z1(1)
            z0(2:end+1,1)=z0(1:end);z0(1)=0; 
        end
        if length(z0)>length(z1), z0(end)=[]; end
        l=sort(z1-z0); 
        l1=l(find(l>lmin));
    end
    DET(w-j+1)=sum(l1)/sum(X);
    L(w-j+1)=mean(l1);
    RR(w-j+1)=sum(X)/length(X);

    j=j+1;

end


L(find(isnan(L)))=0;
RR(find(isnan(RR)))=0;
DET(find(isnan(DET)))=0;

if nargout, XCF=xcf(x1,x2,w,1); end

if ~nargout
subplot(2,2,1)
clim=1;
xcf(x,y,w)
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
Lmax = max(L); if ~Lmax, Lmax = 1; end
axis([-w w 0 Lmax]),ylabel('Averaged Line Length'),grid on
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







% compute the diagonalwise RQA
clear DET L RR, j = 0;
for i = 0:w,waitbar((i+w)/(2*w))
clear z z0 z1

    if m > 1

        switch(method)

        %%%%%%%%%%%%%%%%% local CRP, fixed distance

          case 1
        %%%%%%%%%%%%%%%%% maximum norm
            s = max(abs(-x2(1:NX-i,:) - y2(i+1:NY,:))');
          case 2
        %%%%%%%%%%%%%%%%% euclidean norm
            errcode=112;
            s = sqrt(sum((-x2(1:NX-i,:) - y2(i+1:NY,:)).^2, 2));
          case 3
        %%%%%%%%%%%%%%%%% minimum norm
            errcode=113;
            s = sum(abs(-x2(1:NX-i,:) - y2(i+1:NY,:))');
        end
    else
        s = abs(-x2(1:NX-i) - y2(i+1:NY));
    end
      
    X = s(:) < e;

    if sum(X) == length(X), l1=length(X); else
        z=diff(X); z0 = [];
        if ~isempty(find(~(z-1))),z0(:,1)=find(~(z-1));else,z0(1:length(X))=0;end, 
        if ~isempty(find(~(z+1))),z1=find(~(z+1));else,z1(1:length(X))=0;end

        if z0(1)>z1(1)
            z0(2:end+1,1)=z0(1:end);z0(1)=0; 
        end
        if length(z0)>length(z1), z0(end)=[]; end
        l=sort(z1-z0); l1=l(find(l>lmin));
    end
    DET(j+w+1)=sum(l1)/sum(X);
    L(j+w+1)=mean(l1);
    RR(j+w+1)=sum(X)/length(X);

    if m > 1

        switch(method)

        %%%%%%%%%%%%%%%%% local CRP, fixed distance

          case 1
        %%%%%%%%%%%%%%%%% maximum norm
            s = max(abs(-x2(i+1:NX,:) - y2(1:NY-i,:))');
          case 2
        %%%%%%%%%%%%%%%%% euclidean norm
            errcode=112;
            s = sqrt(sum((-x2(i+1:NX,:) - y2(1:NY-i,:)).^2, 2));
          case 3
        %%%%%%%%%%%%%%%%% minimum norm
            errcode=113;
            s = sum(abs(-x2(i+1:NX,:) - y2(1:NY-i,:))');
        end
    else
        s = abs(-x2(i+1:NX,:) - y2(1:NY-i,:));
    end
    
    X = s(:) < e;

    if sum(X) == length(X), l1=length(X);else
        z=diff(X); z0 = [];
        if ~isempty(find(~(z-1))),z0(:,1)=find(~(z-1));else,z0(1:length(X))=0;end, 
        if ~isempty(find(~(z+1))),z1=find(~(z+1));else,z1(1:length(X))=0;end

        if z0(1)>z1(1)
            z0(2:end+1,1)=z0(1:end);z0(1)=0; 
        end
        z1 = z1(:); z0 = z0(:);

        if length(z0)>length(z1), z0(end)=[]; end
        l=sort(z1-z0); l1=l(find(l>lmin));
    end
    
    DET(w-j+1)=sum(l1)/sum(X);
    L(w-j+1)=mean(l1);
    RR(w-j+1)=sum(X)/length(X);

    j=j+1;

end
if ishandle(hw), delete(hw), end

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

