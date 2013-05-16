function xout=crp_big(varargin)
%CRP_BIG   Creates a cross recurrence plot/ recurrence plot.
%    CRP_BIG(X [,Y] [,param1,param2,...]) creates a cross 
%    recurrence plot/ recurrence plot. Results can be 
%    stored into the workspace. In contrast to CRP, long 
%    data series can be used. Results can be stored into 
%    the workspace.
%
%    R=CRP_BIG(X,M,T,E) uses the dimension M, delay T 
%    and the size of neighbourhood E and creates a recurrence 
%    plot of X.
%    
%    R=CRP_BIG(X,Y,'order') creates an order matrix plot
%    with normalization of the data.
%    
%    R=CRP_BIG(X,Y,'distance','nonormalize') creates a 
%    distance coded matrix plot without normalization
%    of the data.
%    
%    Allows to change the parameters interactively by 
%    using a GUI.
%
%    The source-data X and test-data Y can be one- or 
%    a two-coloumn vectors (then, in the first column 
%    have to be the time); if the test-data Y is not
%    specified, a simple recurrence plot is created.
%
%    Parameters: dimension M, delay T and the size of
%    neighbourhood E are the first three numbers after
%    the data series; further parameters can be used
%    to switch between various methods of finding the
%    neighbours of the phasespace trajectory, to suppress
%    the normalization of the data and to suppress the 
%    GUI (useful in order to use this programme by other 
%    programmes).
%
%    Methods of finding the neighbours/ of plot.
%      maxnorm     - Maximum norm.
%      euclidean   - Euclidean norm.
%      minnorm     - Minimum norm.
%      nrmnorm     - Euclidean norm between normalized vectors
%                    (all vectors have the length one).
%      maxnorm     - Maximum norm, fixed recurrence rate.
%      fan         - Fixed amount of nearest neighbours.
%      inter       - Interdependent neighbours.
%      omatrix     - Order matrix.
%      opattern    - Order patterns recurrence plot.
%      distance    - Distance coded matrix (global CRP, Euclidean norm).
%
%    Normalization of the data series.
%      normalize   - Normalization of the data.
%      nonormalize - No normalization of the data.
%
%    Suppressing the GUI.
%      gui         - Creates the GUI and the output plot.
%      nogui       - Suppresses the GUI and the output plot.
%      silent      - Suppresses all output.
%
%    Parameters not needed to be specified.
%
%    Current limitation: for higher speed in
%    output the whole matrix of the recurrence
%    plot is in the work space - this limits
%    the application of long data series.
%
%    Examples: a = sqrt(100 ^ 2 - (-71:71) .^ 2); b = 1:100;
%              b(101:100 + length(a)) = -(a) + 170;
%              b(end+1:end+100) = 100:-1:1;
%              crp_big(b,1,1,.1,'max')
%
%    See also CRP, CRP2, JRP, TAUCRP and CRQA.

% Copyright (c) 1998-2007 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2008/07/02 11:59:43 $
% $Revision: 5.12 $
%
% $Log: crp_big.m,v $
% Revision 5.12  2008/07/02 11:59:43  marwan
% bug fix for logical data vectors
%
% Revision 5.11  2007/07/18 17:18:44  marwan
% integer values in the arguments supported
%
% Revision 5.10  2007/05/15 17:33:13  marwan
% new neighbourhood criterion: fixed RR
%
% Revision 5.9  2007/03/29 13:17:51  marwan
% uint8 of order patterns and order matrix
%
% Revision 5.8  2006/10/24 14:16:16  marwan
% minor change: sigma in title line of RP shown only for normalised data
%
% Revision 5.7  2006/03/29 13:07:55  marwan
% problems regarding OPRPs and embedding resolved
%
% Revision 5.6  2006/02/06 13:46:17  marwan
% plugin for order patterns recurrence plots supported
%
% Revision 5.5  2005/11/23 07:30:30  marwan
% modified interdependent algorithm
% bug in showing RP fixed
%
% Revision 5.4  2005/11/09 08:58:13  marwan
% bug fix in interdependent neighbours method
%
% Revision 5.3  2005/06/15 15:16:00  marwan
% RP matrix transposed for "interdependent neighbourhood"
%
% Revision 5.2  2005/04/15 09:02:32  marwan
% minor bugfix in plugin section
%
% Revision 5.1  2005/04/08 09:52:27  marwan
% plugin added
%
% Revision 4.8.1.6  2005/03/16 12:21:54  marwan
% add hint in help text for joint recurrence plots
%
% Revision 4.8.1.5  2005/03/16 11:19:02  marwan
% help text modified
%
% Revision 4.8.1.4  2005/01/25 10:14:01  marwan
% bug in order pattern RP and order matrix (matrix transposed)
%
% Revision 4.8.1.3  2004/12/23 07:49:03  marwan
% bug in order patterns RP fixed (empty order patterns)
%
% Revision 4.8.1.2  2004/11/15 12:52:35  marwan
% bug fix in choice of neighbourhood (enabling and disabling parameters)
%
% Revision 4.8.1.1  2004/11/12 08:39:03  marwan
% bug fix in order patterns representation
%
% Revision 4.8  2004/11/11 12:18:13  marwan
% order patterns recurrence plot added
%
% Revision 4.7  2004/11/10 07:04:40  marwan
% initial import
%
%
% This program is part of the new generation XXII series.
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.




warning off
global errcode props nonorm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% programme properties

errcode=0;
init_properties
hCRP=[];hCtrl=[];nogui=[];obj=[];mflag=[];
set(0,'ShowHidden','On')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check and read the input
error(nargchk(1,9,nargin));
if nargout>1, error('Too many output arguments'), end

check_meth={'ma','eu','mi','nr','rr','fa','in','om','op','di'}; 	% maxnorm, euclidean, nrmnorm,  fan, distance
check_norm={'non','nor'};				% nonormalize, normalize
check_gui={'gui','nog','sil'};				% gui, nogui, silent

% transform any int to double
intclasses = {'uint8';'uint16';'uint32';'uint64';'int8';'int16';'int32';'int64';'logical'};
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

if isnumeric(varargin{1})==1 		% read commandline input
   varargin{9}=[];
   i_double=find(cellfun('isclass',varargin,'double'));
   i_char=find(cellfun('isclass',varargin,'char'));

   % check the text input parameters for method, gui and normalization
   temp_meth=0;
   temp_norm=0;
   temp_gui=0;
   if ~isempty(i_char)
      for i=1:length(i_char), 
         varargin{i_char(i)}(4)='0';
         temp_meth=temp_meth+strcmpi(varargin{i_char(i)}(1:2),check_meth'); 
         temp_norm=temp_norm+strcmpi(varargin{i_char(i)}(1:3),check_norm'); 
         temp_gui=temp_gui+strcmpi(varargin{i_char(i)}(1:3),check_gui'); 
      end
      method=min(find(temp_meth));
      nonorm=min(find(temp_norm))-1;
      nogui=min(find(temp_gui))-1;
      if isempty(method), method=1; end
      if isempty(nonorm), nonorm=1; end
      if isempty(nogui), nogui=0; end
      if method>length(check_meth), method=length(check_meth); end
      if nonorm>1, nonorm=1; end
      if nogui>2, nogui=2; end
   else
      method=1; nonorm=1; nogui=0;
   end
   if nogui==0 & nargout>0, nogui=1; end

   % get the parameters for creating RP
     if max(size(varargin{1}))<=3
        error('To less values in data X.')
     end
     x=double(varargin{1});
%     if isempty(varargin{2}) | ~isnumeric(varargin{2}), y=x; end
%     if ~isempty(varargin{2}) & isnumeric(varargin{2}), y=double(varargin{2}); end
     if isempty(varargin{2}) | ~isnumeric(varargin{2}), y=x; 
     else, y=double(varargin{2}); end

     if (isnumeric(varargin{2}) & max(size(varargin{2}))==1) | ~isnumeric(varargin{2})
       y=x;
       if ~isempty(varargin{i_double(2)}), m=varargin{i_double(2)}(1); else m=1; end
       if ~isempty(varargin{i_double(3)}), t=varargin{i_double(3)}(1); else t=1; end
       if ~isempty(varargin{i_double(4)}), e=varargin{i_double(4)}(1); else e=.1; end
     else
       if ~isempty(varargin{i_double(3)}), m=varargin{i_double(3)}(1); else m=1; end
       if ~isempty(varargin{i_double(4)}), t=varargin{i_double(4)}(1); else t=1; end
       if ~isempty(varargin{i_double(5)}), e=varargin{i_double(5)}(1); else e=.1; end
     end
     if e<0, e=1; disp('Warning: The threshold size E cannot be negative and is now set to 1.'), end
     if t<1, t=1; disp('Warning: The delay T cannot be smaller than one and is now set to 1.'), end
     t=round(t); m=round(m); mflag=method;
     if m < 1, m = 1; end
     if t < 1, t = 1; end
     if method==8 & m > 1, 
         m=1; 
         disp('Warning: For order matrix a dimension of one is used.')
     end
     if method==9 & m == 1, 
         m=2; 
         disp(['Warning: For order patterns recurrence plots the dimension must',10,...
              'be larger than one. ',...
              'Embedding dimension is set to ',num2str(m),'.'])
     end
     action='init';

    Nx=length(x); Ny=length(y);
    NX=Nx-t*(m-1);NY=Ny-t*(m-1);  
    if (NX<1 | NY<1) & nogui==0;
         errordlg('The embedding vectors cannot be created. Dimension M and/ or delay T are to big. Please use smaller values.','Dimension/ delay to big')
         waitforbuttonpress
    end
  if size(x,1)<size(x,2), x=x'; end
  if size(y,1)<size(y,2), y=y'; end

  if ~isempty(find(isnan(x)))
     disp('Warning: NaN detected (in first variable) - will be cleared.')
     for k=1:size(x,2),  x(find(isnan(x(:,k))),:)=[]; end
  end
  if ~isempty(find(isnan(y)))
     disp('Warning: NaN detected (in second variable) - will be cleared.')
     for k=1:size(y,2),  y(find(isnan(y(:,k))),:)=[]; end
  end

  if size(x,2)>=2
     xscale=x(:,1); 
     if ~isempty(find(diff(xscale)<0)), error('First column of the first vector must be monotonically non-decreasing.'),end
     if nonorm==1, x=(x(:,2)-mean(x(:,2)))/std(x(:,2)); else x=x(:,2); end
  else
     if nonorm==1, x=(x-mean(x))/std(x); end
     xscale=(1:length(x))'; 
  end
  if size(y,2)>=2
     yscale=y(:,1); 
     if ~isempty(find(diff(yscale)<0)), error('First column of the second vector must be monotonically non-decreasing.'),end
     if nonorm==1, y=(y(:,2)-mean(y(:,2)))/std(y(:,2)); else y=y(:,2); end
  else
      if nonorm==1, y=(y-mean(y))/std(y); end
      yscale=(1:length(y))';
  end
  ds=eye(m);


else 			%  read input from the GUI
  action=varargin{1};
  nogui=0;
  h=get(gcf,'Name');h=h(findstr(h,'(')+1:findstr(h,')')-1);
  hCRP=findobj('Name',['Cross Recurrence Plot (' h ')']);
  hCtrl=findobj('Name',['Control (' h ')']);
  h=str2num(h);
  xshuttle=get(findobj('Parent',hCRP,'Tag','DataPlot2'),'UserData');
  if ~isempty(xshuttle)
    xscale=xshuttle(:,1);
    x=xshuttle(:,2);
    yshuttle=get(findobj('Parent',hCRP,'Tag','DataPlot1'),'UserData');
    yscale=yshuttle(:,1);
    y=yshuttle(:,2);

    if ~isempty(hCtrl)
      if get(findobj('Tag','Unthresh','Parent',hCtrl),'Value')
         mflag=length(check_meth);
      else   
         mflag=get(findobj('Tag','Method','Parent',hCtrl),'Value');
      end
      m=get(findobj('Tag','Dim','Parent',hCtrl),'Value');
      t=str2num(get(findobj('Tag','Delay','Parent',hCtrl),'String'));
      e=str2num(get(findobj('Tag','Size','Parent',hCtrl),'String'));
      if e<0, e=1; 
          errordlg('The threshold size E cannot be negative.','Threshold size to small')
          waitforbuttonpress
  	  set(findobj('Tag','Size','Parent',hCtrl),'String','1')
          action='';
      end
      if t<1, t=1; 
          errordlg('The delay T cannot be smaller than one.','Delay to small')
          waitforbuttonpress
  	  set(findobj('Tag','Delay','Parent',hCtrl),'String','1')
          action='';
      end
      ds = get(findobj('Tag','Dim','Parent',hCtrl),'UserData');
      if mflag==8 | mflag==9 | mflag==length(check_meth)
         set(findobj('Tag','Size','Parent',hCtrl),'enable','off');
         set(findobj('Tag','Sizetext','Parent',hCtrl),'enable','off');
      else
         set(findobj('Tag','Size','Parent',hCtrl),'enable','on');
         set(findobj('Tag','Sizetext','Parent',hCtrl),'enable','on');
      end
      nonorm = get(findobj('Tag','nonorm','Parent',hCtrl),'UserData');
      
      if mflag==8 & m > 1, 
          m=1; ds = 1;
          disp('Warning: For order matrix a dimension of one is used.')
      end
      if mflag==9 & m == 1, 
          m=2;
          errordlg(['For order patterns recurrence plots the',10,'dimension must be larger than one.'],'Dimension too small')
          set(findobj('Tag','Dim','Parent',hCtrl),'Value',m)
      end
      Nx=length(x); Ny=length(y);
      NX=Nx-t*(m-1);NY=Ny-t*(m-1);  
      if (NX<1 | NY<1) & strcmpi(action,'apply');
         errordlg('The embedding vectors cannot be created. Dimension M and/ or delay T are to big. Please use smaller values.','Dimension/ delay to big')
         waitforbuttonpress
         action='';
      end
    end
  end
  clear xshuttle yshuttle
  cm_old=get(hCRP,'Colormap');
  cm=[{hsv(256)}; {hot(256)}; {gray(256)};...
    {french(256)}; {bone(256)}; {copper(256)}; {pink(256)};...
    {flag(256)}; {lines(256)}; {colorcube(256)};...
    {jet(256)};  {prism(256)}; {cool(256)};...
    {autumn(256)}; {spring(256)}; {winter(256)};...
    {summer(256)}; {flipud(gray(2))}; {flipud(cm_old)}];

  if isempty(findobj('Tag','CRPFig')) & nogui==0
    action='init';
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% splash the GPL

splash_gpl('crp')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% nogui

if nogui>0
   hCRP=9999;
   if nogui~=2 
       tx(1)={'maximum norm'}; 
       tx(2)={'Euclidean norm'}; 
       tx(3)={'minimum norm'}; 
       tx(4)={'Euclidean norm of normalized distance'}; 
       tx(5)={'maximum norm fixed RR'};
       tx(6)={'fixed amount of nearest neighbours'};
       tx(7)={'interdependent neighbours'};
       tx(8)={'order matrix'};
       tx(9)={'order pattern'};
       tx(10)={'distance plot'};
       disp(['use method: ', char(tx(method))]);
       if nonorm==1, disp('normalize data'); else disp('do not normalize data'); end
   end
   action='compute';
   if (NX<1 | NY<1)
         disp('Warning: The embedding vectors cannot be created.')
	 disp('Dimension M and/ or delay T are to big. Please use smaller values.')
         action='';
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% switch routines


try
switch(action)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initialization

  case 'init'

  errcode=1;
  xshuttle(:,1)=xscale;
  yshuttle(:,1)=yscale;
  xshuttle(:,2)=x;
  yshuttle(:,2)=y;
  ds=eye(m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% create GUI

  errcode=2;
  root_ud=get(0,'UserData'); 
  if isstruct(root_ud)
    if isfield(root_ud,'crp')
      if ~isempty(root_ud.crp)
        root_ud.crp=[root_ud.crp max(root_ud.crp)+1];
      else
        root_ud.crp=1;
      end
      h=num2str(root_ud.crp(end));
    else
      root_ud.crp=1;
      h=num2str(1);
    end
  else
    root_ud.old=root_ud;
    root_ud.crp=1;
    h=num2str(1);
  end
  set(0,'UserData',root_ud)
    
  %%%%%%%%%%%%%%%%% CRP Figure

   h_axes=create_CRPfig(h,xshuttle,yshuttle);
   if nargout>0, xout=haxes; end	    

  %%%%%%%%%%%%%%%%% Control Figure

  errcode=3;
  create_Ctrlfig('crp_big',h,ds,m,t,e,method)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  vectorswitch

  case 'vectorswitch'
  vectorswitch
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  fit dimension display

  case 'fitdim'

  errcode=5;
	ds2=eye(m); 
	if m>length(ds)
	  ds2(1:length(ds),1:length(ds))=ds;
	else
	  ds2=ds(1:m,1:m);
	end
	ds=ds2; clear ds2
	set(findobj('Tag','Dim','Parent',gcf),'UserData', ds);
	for i=1:20
	   if i>m, set(findobj('Tag',['DimShift' num2str(i)],'Parent',gcf), 'Enable', 'off');
	     else, set(findobj('Tag',['DimShift' num2str(i)],'Parent',gcf), 'Enable', 'on'); end
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% unthresh

  case 'unthresh'

  errcode=6;
  switch_unthresholded(hCRP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stretch

  case 'stretch'

  errcode=7;
  stretch(hCRP,xscale,yscale,Nx,Ny)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change colormap

  case 'log'
 
  errcode=82;
  change_colormapscale(hCRP,cm)


  case 'colormap'
 
  errcode=81;
  change_colormap(hCtrl,hCRP,cm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% store

  case 'store'
  
  errcode=9;
  X=get(findobj('Tag','CRPData','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'UserData');
  if isempty(X)
     warndlg('The CRP matrix is still empty. Please start the computation before storing.','Empty CRP')
     waitforbuttonpress
  else
     assignin('base','X', double(X))
  end
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  handlevisON

  case 'handlevisON'
  
  set(hCRP, 'HandleVis','on')
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% close

  case 'close'
  errcode=101;
  close_all  


  case 'smartclose'
  errcode=102;
  if ishandle(hCRP) & ishandle(hCtrl)
    smart_close(hCRP,hCtrl)
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% computation
  
  case 'compute'

  errcode=11;
  txt_cross='Cross ';
  if size(x)==size(y) 
    if x==y, txt_cross=''; end
  end
  if nogui==0
    setptr([hCRP,hCtrl],'watch')
    h=findobj('Tag','Apply','Parent',hCtrl);
    obj.children=get(hCtrl,'Children');
    obj.enable=get(obj.children,'Enable'); 
    set(obj.children,'Enable','off')
    set(h(1),'String','Stop',...
             'ToolTip','Stops the computation.',...
             'Enable','on',...
	     'Callback','set(gcbo,''String'',''Stopped'')')
  end

  % check if plugin exist and is executable
  [plugin_exist, plugin_name, plugin_path] = is_crp_plugin;
  if plugin_exist & ( mflag < 4 | mflag == 9 | mflag == 10 ) & length(x) == length(y) & ~ispc % if plugin exist and method is MAX, MIN, EUC ord DIS
      if nogui == 1, disp('(plugin used)'), end
  end
  

  maxN=sqrt(m)*max(NX,NY);
  if maxN>800
    n=round(log(maxN/350)/log(2)); % orig war 150
    if plugin_exist & ( mflag < 4 | mflag == 8 | mflag == 9 ) & length(x) == length(y), n=round(log(maxN/550)/log(2)); end
    maxN=maxN/(2^n);
    maxN=round(maxN);NX0=NX;NY0=NY;smax=sqrt(m)*(max(abs(x))+max(abs(y)));
  else
    maxN=round(sqrt(m)*800);NX0=NX;NY0=NY;smax=sqrt(m)*(max(abs(x))+max(abs(y)));
  end
  s_norm=sqrt(m)*abs(max(x(:))-min(x(:)));
  set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'Visible','on')
  set(findobj('Tag','CRPData','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'Visible','off')
  set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Create CRP matrix'),drawnow

  X=uint8(zeros(NY+(m-1)*t,NX+(m-1)*t));
  set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Build Embedding Vectors'),drawnow
    if m*Nx>maxN
      ax=ceil(m*Nx/maxN);
      Nx2=ceil(Nx/ax);
      x(Nx+1:Nx2*ax+t*(m-1))=x(Nx);
    else
      ax=1; Nx2=Nx;
    end
    x(Nx+1:Nx2*ax+t*(m-1))=x(Nx);
    x0=x;

    if m*Ny>maxN
      ay=ceil(m*Ny/maxN);
      Ny2=ceil(Ny/ay);
    else
      ay=1; Ny2=Ny;
    end
    y(Ny+1:Ny2*ay+t*(m-1))=y(Ny);
    y0=y;
   
       
  if ~nonorm, unit = ''; else unit = '\sigma'; end

  set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Compute CRP'),drawnow
  if nogui==0
%    set(hCRP,'Renderer','OpenGL')
    h1=axes('Parent',hCRP,'pos',[.3 .34 .4 .02],'xtick',[],'ytick',[],'box','on','NextPlot','add');
    h2=fill([0 0 0 0],[0 0 1 1],[.2 .3 1],'Parent',h1,'LineStyle','None');
    set(h1,'xtick',[],'ytick',[],'box','on','xlim',[0 1],'layer','top')
  end
  k=0;
  for i=1:ax,
        if nogui==0
          h=findobj('Tag','Apply','Parent',hCtrl);
          if strcmpi(get(h(1),'String'),'stopped')
            set(h(1),'String','Apply',...
     	             'ToolTip','Starts the computation - be patient.',...
	             'Callback','crp_big compute')
            for i=1:length(obj.enable), set(obj.children(i),'Enable',obj.enable{i}); end
            break
          end
        end

        % make substring of data
        x=x0(1+Nx2*(i-1):Nx2+Nx2*(i-1)+t*(m-1));
          for j=1:ay,
                y=y0(1+Ny2*(j-1):Ny2+Ny2*(j-1)+t*(m-1));
	            k=k+1;

                Nx=length(x); Ny=length(y);
                NX=Nx-t*(m-1);NY=Ny-t*(m-1);  
                if nogui==0
                  if k/10==fix(k/10), set(h2,'XData',[0 k/(ax*ay) k/(ax*ay) 0]), drawnow, end
                end
                
                % embedding vectors
                ii=(1:NX)';jj=0:t:0+(m-1)*t;
                ii=reshape(ii(:,ones(1,m))+jj(ones(NX,1),:),m*NX,1);
                x1=x(ii);
                x2=reshape(x1,NX,m);

                ii=(1:NY)';jj=0:t:0+(m-1)*t;
                ii=reshape(ii(:,ones(1,m))+jj(ones(NY,1),:),m*NY,1);
                y1=y(ii);
                y2=reshape(y1,NY,m);

                y2=y2*ds; % switch vectors

                %%

                if plugin_exist & ( mflag < 4 | mflag == 8 | mflag == 9 ) & length(x) == length(y) % if plugin exist and method is MAX, MIN, EUC ord DIS

                    [X0 matext] = crp_plugin(x2, y2, 1, 1, e, mflag, hCRP, plugin_path, 1);
                    X(1+Ny2*(j-1):Ny2+Ny2*(j-1),1+Nx2*(i-1):Nx2+Nx2*(i-1))=X0;
                    X(NY0+1:end,:)=[];
                    X(:,NX0+1:end)=[];
                    NX=NX0; NY=NY0;

                else
                % use builtin implementation


                  [NX, mx] = size(x2);
                  [NY, my] = size(y2);

                  clear jx jy

                  switch(mflag)

                      %%%%%%%%%%%%%%%%% local CRP, fixed distance

                      case {1,2,3,5}

                      errcode=111;
                      px = permute(x2, [ 1 3 2 ]);
                      py = permute(y2, [ 3 1 2 ]);
                      s1 = px(:,ones(1,NY),:) - py(ones(1,NX),:,:);
                      switch mflag
                        case 1
                      %%%%%%%%%%%%%%%%% maximum norm
                          s = max(abs(s1),[],3);
                          matext=[num2str(round(100*e)/100) unit ' (fixed distance maximum norm)'];
                        case 5
                      %%%%%%%%%%%%%%%%% maximum norm, fixed RR
                          errcode=115;
                          s = max(abs(s1),[],3);
                          ss = sort(s(:));
                          idx = ceil(e * length(ss));
                          e = ss(idx);
                          matext=[num2str(round(100*e)/100) '\sigma (fixed distance maximum norm, fixed RR)'];
                        case 2
                      %%%%%%%%%%%%%%%%% euclidean norm
                          errcode=112;
                          s = sqrt(sum(s1.^2, 3));
                          matext=[num2str(round(100*e)/100) unit ' (fixed distance euclidean norm)'];
                        case 3
                      %%%%%%%%%%%%%%%%% minimum norm
                          errcode=113;
                          s = sum(abs(s1), 3);
                          matext=[num2str(round(100*e)/100) unit ' (fixed distance minimum norm)'];
                      end

                      X1=255*s/max(s(:))<(255*e/max(s(:)));
                      X0=(uint8(X1))'; clear s s1 x1 y1 px py X1
                      X(1+Ny2*(j-1):Ny2+Ny2*(j-1),1+Nx2*(i-1):Nx2+Nx2*(i-1))=X0;
                      X(NY0+1:end,:)=[];
                      X(:,NX0+1:end)=[];
                      NX=NX0; NY=NY0;
                    %  matext=[num2str(round(100*e)/100) '\sigma (fixed distance maximum norm)'];



                      %%%%%%%%%%%%%%%%% local CRP, normalized distance euclidean norm

                      case 4

                      errcode=114;
                    %  set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String',['Normalize Embedding Vectors ',num2str(ax*(i-1)+j)]),drawnow
                      Dx=sqrt(sum(((x2(:,:)).^2)'))';
                      Dy=sqrt(sum(((y2(:,:)).^2)'))';
                      x1=x2./repmat(Dx,1,m);x2=x1;
                      y1=y2./repmat(Dy,1,m);y2=y1; clear Dx Dy y1 x1


                      px = permute(x2, [ 1 3 2 ]);
                      py = permute(y2, [ 3 1 2 ]);
                      s1 = px(:,ones(1,NY),:) - py(ones(1,NX),:,:);
                      s = sqrt(sum(s1.^2, 3));

                      X0=(uint8(255*s/max(s(:)))<(255*e/max(s(:))))'; clear s s1 x1 y1 px py

                      X(1+Ny2*(j-1):Ny2+Ny2*(j-1),1+Nx2*(i-1):Nx2+Nx2*(i-1))=X0;
                      X(NY0+1:end,:)=[];
                      X(:,NX0+1:end)=[];
                      NX=NX0; NY=NY0;
                      matext=[num2str(round(100*e)/100) unit ' (normalized distance euclidean norm)'];


                      %%%%%%%%%%%%%%%%% local CRP, fixed neigbours amount

                      case 6

                      errcode=116;
                      if e>=1 
                        e=round(e)/100;
                        txt=['The value for fixed neigbours amount has to be smaller '...
                             'than one. Continue the computation with a value of ' ...
	                     num2str(e)];
                        if nogui==0
                          warndlg(txt,'Threshold value mismatch');
                          drawnow
                          waitforbuttonpress
                          set(findobj('Tag','Size','Parent',gcf),'String',num2str(e))
                        end
                      end


                      px = permute(x2, [ 1 3 2 ]);
                      py = permute(y2, [ 3 1 2 ]);
                      s1 = px(:,ones(1,NY),:) - py(ones(1,NX),:,:);
                      s = sqrt(sum(s1.^2, 3));

                      mine=round(NY*e);
                      [SS, JJ]=sort(s'); JJ=JJ';

                      X1(NX*NY)=uint8(0); X1(JJ(:,1:mine)+repmat([0:NY:NX*NY-1]',1,mine))=uint8(1);

                      X0=reshape(X1,NY,NX); clear X1 SS JJ s px py
                      X(1+Ny2*(j-1):Ny2+Ny2*(j-1),1+Nx2*(i-1):Nx2+Nx2*(i-1))=X0;
                      X(NY0+1:end,:)=[];
                      X(:,NX0+1:end)=[];
                      matext=[num2str(round(1000*e)/10) '% (fixed neighbours amount)'];
                      NX=NX0; NY=NY0;


                      %%%%%%%%%%%%%%%%% local CRP, interdependent neigbours 

                      case 7

                      errcode=117;
                      if e>=1 
                        e=round(e)/100;
                        txt=['The value for fixed neigbours amount has to be smaller '...
                             'than one. Continue the computation with a value of ' ...
	                     num2str(e)];
                        if nogui==0
                          warndlg(txt,'Threshold value mismatch');
                          drawnow
                          waitforbuttonpress
                          set(findobj('Tag','Size','Parent',gcf),'String',num2str(e))
                        end
                      end

                      px = permute(x2, [ 1 3 2 ]);
                      py = permute(y2, [ 1 3 2 ]);
                      px2 = permute(x2, [ 3 1 2 ]);
                      py2 = permute(y2, [ 3 1 2 ]);

                      s1 = px(:,ones(1,NX),:) - px2(ones(1,NX),:,:);
                      sx = sqrt(sum(s1.^2, 3));
                      s1 = py(:,ones(1,NY),:) - py2(ones(1,NY),:,:);
                      sy = sqrt(sum(s1.^2, 3));

                      mine=round(min(NX,NY)*e);
                      [SSx, JJx]=sort(sx); %SSx(1,:)=[]; JJx(1,:)=[];
                      [SSy, JJy]=sort(sy); %SSy(1,:)=[]; JJy(1,:)=[];
                      eey=mean(SSy(mine:mine+1,:));
                      eex=mean(SSx(mine:mine+1,:));
                      X0=uint8(zeros(NY,NX));
                      for ii=1:min(NX,NY),
                        JJx((JJx(1:mine,ii)>min(NX,NY)),ii)=i;
                        X0(ii,JJx(1:mine,ii)) = sy(ii,JJx(1:mine,ii))<=eey(ii);
                        JJy((JJy(1:mine,ii)>min(NX,NY)),ii)=i;
                        X0(ii,JJy(1:mine,ii)) = sx(ii,JJy(1:mine,ii))<=eex(ii);
                      end

                      X(1+Ny2*(j-1):Ny2+Ny2*(j-1),1+Nx2*(i-1):Nx2+Nx2*(i-1))=X0';
                      X(NY0+1:end,:)=[];
                      X(:,NX0+1:end)=[];
                      NX=NX0; NY=NY0;
                      clear X1 SS* JJ* s sx sy s1 px py eex eey
                      X = X';
                      matext=[num2str(round(1000*e)/10) '% (interdependent neighbours)'];


                      %%%%%%%%%%%%%%%%% order matrix

                      case 8

                      errcode=118;

                      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Compute Order Matrix'),drawnow
                      px = permute(x2, [ 1 3 2 ]);
                      py = permute(y2, [ 3 1 2 ]);
                      X0 = uint8(px(:,ones(1,NY),:) >= py(ones(1,NX),:,:) - e);
                      if check_stop(hCRP,hCtrl,nogui,obj), return, end
                      X(1+Ny2*(j-1):Ny2+Ny2*(j-1),1+Nx2*(i-1):Nx2+Nx2*(i-1))=X0';
                      X(NY0+1:end,:)=[];
                      X(:,NX0+1:end)=[];
                      NX=NX0; NY=NY0;

                      matext='';

                      %%%%%%%%%%%%%%%%% order pattern recurrence plot

                      case 9

                      errcode=119;

                      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Compute Order Patterns Matrix'),drawnow

                      % create a order pattern test
                      cmdStr = '';
                      for i3=2:m;
                          cmdStr = [cmdStr, ' permX(:,', num2str(i3-1) ,') < permX(:,', num2str(i3), ') + eps'];
                          if i3 < m, cmdStr = [cmdStr, ' &']; end
                      end
                      if m==1
                          cmdStr = '1'; 
                          disp('Warning: No order patterns for dimension one; please use higher dimension!')
                      end


                      % order patterns by permutation of the set of values
                      clear patt*
                      for i3=1:length(x2)
                          permX=perms(x2(i3,:));
                          orderPattern = find(eval(cmdStr));
                          if isempty(orderPattern) orderPattern = 0; end
                          pattX(i3) = orderPattern(1);
                      end
                      for i3=1:length(y2)
                          permX=perms(y2(i3,:));
                          orderPattern = find(eval(cmdStr));
                          if isempty(orderPattern) orderPattern = 0; end
                          pattY(i3) = orderPattern(1);
                      end

                      px = permute(pattX', [ 1 3 2 ]);
                      py = permute(pattY', [ 3 1 2 ]);
                      X0 = uint8(px(:,ones(1,length(y2)),:) == py(ones(1,length(x2)),:,:));

                      X(1+Ny2*(j-1):Ny2+Ny2*(j-1),1+Nx2*(i-1):Nx2+Nx2*(i-1))=X0';

                      X(NY0+1:end,:)=[];
                      X(:,NX0+1:end)=[];
                      NX=NX0; NY=NY0;

                      matext='';


                      %%%%%%%%%%%%%%%%% global CRP

                      case length(check_meth)

                      errcode=110+length(check_meth);
                      px = permute(x2, [ 1 3 2 ]);
                      py = permute(y2, [ 3 1 2 ]);
                      s1 = px(:,ones(1,NY),:) - py(ones(1,NX),:,:);
                      s = sqrt(sum(s1.^2, 3))';
                      X0=uint8(255*s/s_norm);  

                      X(1+Ny2*(j-1):Ny2+Ny2*(j-1),1+Nx2*(i-1):Nx2+Nx2*(i-1))=X0;
                      X(NY0+1:end,:)=[];
                      X(:,NX0+1:end)=[];
                      NX=NX0; NY=NY0;
                      matext='';


                  end % end switch
              end % end plugin


  
  
%%   
        end % end ay
  end % end ax


  if nogui==0
    for i=1:length(obj.enable), set(obj.children(i),'Enable',obj.enable{i}); end
    set(h(1),'String','Apply',...
             'ToolTip','Starts the computation - be patient.',...
             'Callback','crp_big compute')
    setptr([hCRP,hCtrl],'arrow')
    delete([ h2 h1]),  clear h1 h2
  end
  
  %%%%%%%%%%%%%%%%% show CRP 

  if mflag==6
     disp('Warning! Don''t use these results, because they are not correct.')
  end

  if nogui==0
  set(hCRP,'Renderer','painters')
     Shuttle.hCRP=hCRP;
     Shuttle.hCtrl=hCtrl;
     Shuttle.matext=matext;
     Shuttle.xscale=xscale;
     Shuttle.yscale=yscale;
     Shuttle.mflag=mflag;
     Shuttle.m=m;
     Shuttle.t=t;
     Shuttle.cm=cm;
     Shuttle.txt_cross=txt_cross;
     if isempty(X)
         warn_str = ['Uuups! Empty matrix.',10,'I give up ...'];
         if strcmpi(computer,'GLNXA64')
             warn_str = [warn_str,10,'(Maybe non-appropriate plugin version due to different glibc.)']
         end
         warning(warn_str); 
         return
     end
     show_crp(X,Shuttle)
 else
   if nargout==1, xout=X;end
 end
  
end
try, set(0,props.root), end
warning on
set(0,'ShowHidden','Off')

%%%%%%% error handling

%if 0
catch
  if nargin, xout=NaN; end
  try, if nogui==0
    for i=1:length(obj.enable), set(obj.children(i),'Enable',obj.enable{i}); end
    set(h(1),'String','Apply',...
             'ToolTip','Starts the computation - be patient.',...
             'Callback','crp compute')
    setptr([hCRP,hCtrl],'arrow')
  end, end
  z=whos;x=lasterr;y=lastwarn;in=varargin{1};
  print_error('crp_big',z,x,y,in,mflag,action)
  try, set(0,props.root), end
end
