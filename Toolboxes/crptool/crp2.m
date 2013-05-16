function xout=crp2(varargin)
%CRP2   Creates a cross recurrence plot/ recurrence plot and the LOS.
%    CRP2(X [,Y [,param1,param2,...]) creates a cross recurrence  
%    plot/ recurrence plot from the embedding vectors X and Y. 
%    Results can be stored into the workspace. Further it
%    is possible to estimate the line of synchronization (LOS)
%    in order to get the nonparametric time-relationship between
%    the two considered systems.
%
%    R=CRP2(X,M2,T,E) uses the additionally dimension M2, delay T 
%    and the size of neighbourhood E and creates a recurrence 
%    plot of X.
%    
%    R=CRP2(X,Y,'distance','nonormalize') creates a 
%    distance coded matrix plot without normalization
%    of the data.
%
%    Allows to change the parameters interactively by 
%    using a GUI.
%
%    The embedding dimension M is given by the size of the
%    N x M matrix X and Y; if the matrix Y is not specified, 
%    a simple recurrence plot is created.
%
%    Parameters: additionally dimension M2, delay T and the 
%    size of neighbourhood E are the first three numbers 
%    after the data series; further parameters can be used
%    to switch between various methods of finding the
%    neighbours of the phasespace trajectory, to suppress
%    the normalization of the data and to suppress the 
%    GUI (useful in order to use this programme by other 
%    programmes).
%
%    Methods of finding the neighbours.
%      maxnorm     - Maximum norm.
%      euclidean   - Euclidean norm.
%      minnorm     - Minimum norm.
%      nrmnorm     - Euclidean norm between normalized vectors
%                    (all vectors have the length one).
%      rr          - Maximum norm, fixed recurrence rate.
%      fan         - Fixed amount of nearest neighbours.
%      omatrix     - Order matrix (disabled).
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
%    Examples: a = sin((1:1000) * 2 * pi/200);  % pendulum's location vector
%              b = cos((1:1000) * 2 * pi/200);  % pendulum's velocity vector
%              plot(a,b,'.')
%              crp2(a(1:500),b(1:500),'nonorm','euclidean')
%
%              b = sin(.01 * ([1:1000] * 2 * pi/67) .^2);
%              crp2(b(1:500),a(1:700),3,10,.06,'fan')
%
%    See also CRP, CRP_BIG, JRP, TAUCRP and TRACKPLOT.
%
%    References: 
%    Marwan, N., Thiel, M., Nowaczyk, N.: 
%    Cross Recurrence Plot Based Synchronization of Time Series,
%    Nonlin. Proc. Geophys., 9, 2002.

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 1998-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2012/10/22 14:18:33 $
% $Revision: 5.19 $
%
% $Log: crp2.m,v $
% Revision 5.19  2012/10/22 14:18:33  marwan
% bug fix: normalisation of data when data contains Inf
%
% Revision 5.18  2010/06/29 12:47:30  marwan
% some minor bugs in output and test of time series lengths (of x and y)
%
% Revision 5.17  2009/03/24 08:31:17  marwan
% copyright address changed
%
% Revision 5.16  2008/07/02 11:59:22  marwan
% new norms: DTW and Levenshtein
% bug fix for logical data vectors
%
% Revision 5.15  2008/04/29 14:52:01  marwan
% levenshtein and DTW distance added
%
% Revision 5.14  2007/07/18 17:18:44  marwan
% integer values in the arguments supported
%
% Revision 5.13  2007/05/15 17:33:13  marwan
% new neighbourhood criterion: fixed RR
%
% Revision 5.12  2007/05/15 16:00:38  marwan
% minor change in outfit
%
% Revision 5.11  2007/03/29 13:17:51  marwan
% uint8 of order patterns and order matrix
%
% Revision 5.10  2006/10/24 14:16:16  marwan
% minor change: sigma in title line of RP shown only for normalised data
%
% Revision 5.9  2006/03/29 13:07:55  marwan
% problems regarding OPRPs and embedding resolved
%
% Revision 5.8  2006/02/14 11:44:50  marwan
% bug in plugin-call (dim and delay) resolved
%
% Revision 5.7  2006/02/06 15:12:46  marwan
% bug in multi-dimensional embedding solved
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
% Revision 5.3  2005/09/02 08:02:57  marwan
% line fitting algorithm improved (linear interpolation between set points)
%
% Revision 5.2  2005/04/15 09:02:32  marwan
% minor bugfix in plugin section
%
% Revision 5.1  2005/04/08 09:52:27  marwan
% plugin added
%
% Revision 4.8.1.7.2.1  2005/04/08 09:22:13  marwan
% Included plugin
%
% Revision 4.8.1.7  2005/04/06 12:57:50  marwan
% small deviation in threshold application fixed
%
% Revision 4.8.1.6  2005/03/16 12:21:54  marwan
% add hint in help text for joint recurrence plots
%
% Revision 4.8.1.5  2005/03/16 11:19:02  marwan
% help text modified
%
% Revision 4.8.1.4  2004/12/23 07:49:03  marwan
% bug in order patterns RP fixed (empty order patterns)
%
% Revision 4.8.1.3  2004/12/02 13:18:30  marwan
% bug due to missing variable
%
% Revision 4.8.1.2  2004/11/15 12:52:34  marwan
% bug fix in choice of neighbourhood (enabling and disabling parameters)
%
% Revision 4.8.1.1  2004/11/12 08:39:02  marwan
% bug fix in order patterns representation
%
% Revision 4.8  2004/11/11 12:18:12  marwan
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

check_meth={'ma','eu','mi','nr','rr','fa','in','om','op','le','dt','di'}; 	% maxnorm, euclidean, nrmnorm,  fan, distance
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

if isnumeric(varargin{1}) 		% read commandline input
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
     if isempty(varargin{2}) | ~isnumeric(varargin{2}), y=x; else
     y=double(varargin{2}); end

     if (isnumeric(varargin{2}) & max(size(varargin{2}))==1) | ~isnumeric(varargin{2})
       y=x;
       if ~isempty(varargin{i_double(2)}), m0=varargin{i_double(2)}(1); else m0=1; end
       if ~isempty(varargin{i_double(3)}), t=varargin{i_double(3)}(1); else t=1; end
       if ~isempty(varargin{i_double(4)}), e=varargin{i_double(4)}(1); else e=.1; end
     else
       if ~isempty(varargin{i_double(3)}), m0=varargin{i_double(3)}(1); else m0=1; end
       if ~isempty(varargin{i_double(4)}), t=varargin{i_double(4)}(1); else t=1; end
       if ~isempty(varargin{i_double(5)}), e=varargin{i_double(5)}(1); else e=.1; end
     end
     t=round(t); m0=round(m0); mflag=method;
     if e<0, e=1; disp('Warning: The threshold size E cannot be negative and is now set to 1.'), end
     if t<1, t=1; disp('Warning: The delay T cannot be smaller than one and is now set to 1.'), end
     if m0 < 1, m0 = 1; end
     if t < 1, t = 1; end
     if size(x,1)==1, x=x'; end, if size(y,1)==1, y=y'; end 
     m=max([size(x,2) size(y,2)]);

     if method==8 & (m*m0) > 1, 
       m0=1; 
       error(['The neighbourhood criterion ''Oder matrix''',10,'is not implemented - use crp or crp_big instead.'])
     end
     if method==9 & (m*m0) == 1, 
       m0=2; 
         disp(['Warning: For order patterns recurrence plots the dimension must',10,...
              'be larger than one. ',...
              'Embedding dimension is set to ',num2str(m0),'.'])
     end
     action='init';

  if ~isempty(find(isnan(x)))
     disp('NaN detected (in first variable) - will be cleared.')
     for k=1:size(x,2),  x(find(isnan(x(:,k))),:)=[]; end
  end
  if ~isempty(find(isnan(y)))
     disp('NaN detected (in second variable) - will be cleared.')
     for k=1:size(y,2),  y(find(isnan(y(:,k))),:)=[]; end
  end
  if size(x,1) < t*(m-1)+1 | size(y,1) < t*(m-1)+1
     error(['Too less data',10,...
            'Either too much NaN or the number of columns in the vectors do not match.'])
  end

  Nx=size(x,1); Ny=size(y,1);
  NX=Nx-t*(m0-1);NY=Ny-t*(m0-1);
  x0=zeros(Nx,m);y0=zeros(Ny,m);
  x0(1:size(x,1),1:size(x,2))=x; 
  y0(1:size(y,1),1:size(y,2))=y; 

  % normalise the data
  if nonorm == 1,
      for k = 1:size(x0,2)
          idx = find(~isinf(x0(:,k)));
          stdx = std(x0(idx,k));
          meanx = mean(x0(idx,k));
          x(:,k) = (x0(:,k) - meanx) / stdx;
      end
      for k = 1:size(y0,2)
          idy = find(~isinf(x0(:,k)));
          stdy = std(y0(idy,k));
          meany = mean(y0(idy,k));
          y(:,k) = (y0(:,k) - meany) / stdy;
      end
  end

  if ~isempty(find(isnan(x))), for k=1:size(x,2),  x(find(isnan(x(:,k))),:)=[]; end, end
  if ~isempty(find(isnan(y))), for k=1:size(y,2),  y(find(isnan(y(:,k))),:)=[]; end, end

  if size(x,1) < t*(m0-1)+1 | size(y,1) < t*(m0-1)+1
     error(['Too less data',10,...
            'Either too much NaN or the number of columns in the vectors do not match.'])
  end
  if(size(x,2) ~= size(y,2)) 
     error(['Matrix dimensions must agree.'])
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
    x=xshuttle(:,2:end);
    yshuttle=get(findobj('Parent',hCRP,'Tag','DataPlot1'),'UserData');
    y=yshuttle(:,2:end);

    if ~isempty(hCtrl)
      if get(findobj('Tag','Unthresh','Parent',hCtrl),'Value')
         mflag=length(check_meth);
      else   
         mflag=get(findobj('Tag','Method','Parent',hCtrl),'Value');
      end
      m=size(x,2);
      m0=get(findobj('Tag','Dim0','Parent',hCtrl),'Value');
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
      ds = get(findobj('Tag','DimEx','Parent',hCtrl),'UserData');
      if mflag==7 | mflag==8 | mflag==length(check_meth)
         set(findobj('Tag','Size','Parent',hCtrl),'enable','off');
         set(findobj('Tag','Sizetext','Parent',hCtrl),'enable','off');
      else
         set(findobj('Tag','Size','Parent',hCtrl),'enable','on');
         set(findobj('Tag','Sizetext','Parent',hCtrl),'enable','on');
      end
      nonorm = get(findobj('Tag','nonorm','Parent',hCtrl),'UserData');
      
      if mflag==8 & (m * m0) == 1, 
          m0=2;
          errordlg(['For order patterns recurrence plots the',10,'dimension must be larger than one.'],'Dimension too small')
          set(findobj('Tag','Dim0','Parent',hCtrl),'Value',m0)
      end
      Nx=length(x); Ny=length(y);
      NX=Nx-t*(m-1);NY=Ny-t*(m-1);  
      xscale=1:Nx; yscale=1:Ny;
      if (NX<1 | NY<1) & strcmpi(action,'apply');
         errordlg('The embedding vectors cannot be created. Dimension M and/ or delay T are to big. Please use smaller values.','Dimension/ delay to big')
         waitforbuttonpress
         action='';
      end
    end
  end

  temp=get(findobj('Tag','fixp','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'XData');
  if ~isempty(temp)
    for i=1:length(temp), if iscell(temp), fixp(i,1)=temp{i}; else, fixp(i,1)=temp(i); end, end
    temp=get(findobj('Tag','fixp','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'YData');
    for i=1:length(temp), if iscell(temp), fixp(i,2)=temp{i}; else, fixp(i,2)=temp(i); end, end
    fixp = sortrows(fixp,1);
  else 
    fixp=[];
  end
  clear xshuttle yshuttle temp
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

h_meth = findobj('Tag','Method','Parent',hCtrl);
if mflag==7, mflag=2; 
   warndlg(['The neighbourhood criterion ''Interdependent''',10,'is not implemented - use crp or crp_big instead.'],'Neighbourhood'); 
   waitforbuttonpress
   if ~isempty(h_meth)
     set(h_meth,'value',mflag)
   end
   return
end
if mflag==8, mflag=2; 
   warndlg(['The neighbourhood criterion ''Oder matrix''',10,'is not implemented - use crp or crp_big instead.'],'Neighbourhood'); 
   waitforbuttonpress
   if ~isempty(h_meth)
     set(h_meth,'value',mflag)
   end
     set(findobj('Tag','Size','Parent',hCtrl),'enable','on');
     set(findobj('Tag','Sizetext','Parent',hCtrl),'enable','on');
   return
end
method=mflag;

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
       tx(11)={'distance plot'};
       tx(12)={'distance plot'};
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
  xshuttle(:,1)=(1:length(x))';
  xshuttle(:,2:size(x,2)+1)=x;
  yshuttle(:,1)=(1:length(y))';
  yshuttle(:,2:size(y,2)+1)=y;
  ds=eye(m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% create GUI

  errcode=2;
  scr=get(0,'ScreenSize'); 
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

   [h_axes,h_fig]=create_CRPfig(h,xshuttle,yshuttle);
   h0=uicontextmenu('Parent',h_fig);
   set(h_axes,'UIContextMenu',h0)
   h2=line([],[],[],'Parent',h_axes,'visible','off','Tag','Diagonal','color',[1 0 0],...
            'LineWidth',1);
   uimenu(h0, 'Label', 'Set Point', 'Callback', 'crp2 LOSset')

  %%%%%%%%%%%%%%%%% Control Figure

  errcode=3;
  h9=figure('Tag','CRPFig',...			% Control Figure
            'Color',[0.8 0.8 0.8], ...
            'Position',[5*scr(3)/8+10 scr(4)/8 1*scr(3)/6 3*scr(4)/4 ],...
            'NumberTitle','off',...
	        'Name',['Control (' h ')'],...
	        'MenuBar','None',...
	        'DeleteFcn','crp2 handlevisON',...
	        'Resize','Off');
  set(h9,props.window,'Units','Norm')

  h0=uicontrol(props.frame, ... % Frame Embedding
            'Units','Normalized',...
            'UserData',nonorm, ...
            'Tag','nonorm', ...
	  		'Position',[.1 .647 .8 .323]);

  h0=uicontrol(props.text,...	% Text Embedding
            'Units','Normalized', ...
            'FontAngle','italic', ...
            'Position',[.12 .939 .6 .02], ...
            'String','Embedding');
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h1(3) h1(4)])

  h0=uicontrol(props.text,...		% Text Dimension
            'Units','Normalized',...
            'String','Dimension:',...
            'Position',[.16 .905 .35 .02]);
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h1(3) h1(4)])
  if method==8, set(h0,'Enable','Off'); end

  h0=uicontrol(props.text,...		% Text Dimensionvalue
            'Units','Normalized',...
            'Tag','Dim',...
            'String',['x ',num2str(m)],...
            'Position',[.74 .905 .15 .02], ...
            'ToolTip','Embedding dimension.');
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h1(3) h1(4)])
  if method==8, set(h0,'Enable','Off'); end

  h0=uicontrol(props.popup,...		% Input Dimensionfactor
            'Units','Normalized',...
            'Tag','Dim0',...
            'String','1|2|3|4',...
            'Position',[.54 .91 .18 .026], ...
            'Value',m0,...
            'Callback','crp2 dimfit',...
            'ToolTip','Select the embedding dimensionfactor.');
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h2(3) h1(4)])
  if method==8, set(h0,'Enable','Off'); end

  h0=uicontrol(props.text,...		% Text Delay
            'Units','Normalized',...
            'String','Delay:',...
            'Position',[.16 .862 .35 .02]);
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h2(3) h1(4)])
  if method==8, set(h0,'Enable','Off'); end

  h0=uicontrol(props.edit,...		% Input Delay
            'Units','Normalized',...
            'Tag','Delay',...
            'String',t,...
            'Position',[.54 .866 .279 .026],...
            'ToolTip','Insert the embedding delay time.');
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h2(3) h1(4)])
  if m0==1;
    set(h0,'Enable','off')
  end
  if method==8, set(h0,'Enable','Off'); end

  h0=uicontrol(props.text,...		% Text Vector Excluding
            'Units','Normalized',...
            'String','Vector Excluding:',...
            'Position',[.16 .826 .6 .02], ...
            'Tag','DimEx',...
            'UserData', ds);
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h1(3) h1(4)])

  for i=1:3, for j=1:4				% Input Vector Excluding
     h0 = uicontrol(props.button,'Units','normalized', ...
	        'Callback','crp2 exclude', ...
	        'Position',[.18+(j-1)*0.167 0.798-(i-1)*0.034 .14 0.028], ...
	        'String',(i-1)*4+j, ...
	        'ToolTip','Exclude this vector.',...
	        'Tag',['DimShift' num2str((i-1)*4+j)]);
     h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h2(3) h1(4)])
  end, end

  for i=1:12
     if i>m, set(findobj('Tag',['DimShift' num2str(i)]), 'Enable', 'off');				
       else, set(findobj('Tag',['DimShift' num2str(i)]), 'Enable', 'on'); end
  end
  
  h0=uicontrol(props.text,...		% Text Copyright
            'Units','Normalized',...
	        'HorizontalAlignment','center',...
	        'Position',[.15 .65 .7 .06]);
  h1=textwrap(h0,{[char(169),' AGNLD'],'University of Potsdam','1998-2007'});
  set(h0,'String',h1)
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h2(3) h1(4)])

  h0=uicontrol(props.frame, ...  % Frame Neighbourhood
            'Units','Normalized',...
	        'Position',[.1 .449 .8 .18]);

  h0=uicontrol(props.text,...		% Text Neighbourhood
            'Units','Normalized',...
	        'FontAngle','italic', ...
	        'String','Neighbourhood',...
	        'Position',[.12 .598 .6 .02]);
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h1(3) h1(4)])

  h0=uicontrol(props.checkbox, ...		% Button Unthresholded
            'Units','Normalized',...
	        'Position',[.16 .555 .51 .032], ...
	        'String','Unthresholded', ...
	        'CallBack','crp2 unthresh',...
	        'Tag','Unthresh',...
	        'ToolTip','Switch between thresholded and unthresholded CRP.' );

  if method==12, set(h0,'Value',1); end
  if method==8 | method==9, set(h0,'Enable','Off'); end

  h1=uicontrol(props.popup, ...		% Button Unthresholded Scale
            'Units','Normalized',...
	        'Position',[.69 .555 .14 .032], ...
	        'String','1|1/2|1/4|1/6|1/8', ...
	        'Value',1,...
	        'CallBack','crp2 log',...
	        'Enable','off',...
	        'Tag','Log',...
	        'ToolTip','Switch between various scaled CRP.' );

  if method==12, set(h1,'Enable','On'); end

  h2=uicontrol(props.popup,...		% Input Neighbourhood Method
            'Units','Normalized',...
	        'String','Maximum Norm|Euclidean Norm|Minimum Norm|Normalized Norm|Fixed RR|Fixed Amount|Interdependent|Order Matrix|Order Pattern|Levenshtein|DTW',...
	        'Position',[.16 .509 .67 .032],...
	        'CallBack','crp2 unthresh',...
	        'Tag','Method',...
	        'ToolTip','Select the method of finding neighbours.');

  if method==12, set(h2,'Enable','Off'); else, set(h2,'Value',method); end

  h0=uicontrol(props.text,...		% Text Threshold
            'Units','Normalized',...
	        'Tag','Sizetext',...
	        'String','Threshold:',...
	        'Position',[.16 .462 .35 .02]);
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h2(3) h1(4)])

  if method==12 | method==9, set(h0,'Enable','Off'); end


  h0=uicontrol(props.edit,...		% Input Threshold
            'Units','Normalized',...
	        'Tag','Size',...
	        'Position',[.58 .466 .249 .026],...
	        'HorizontalAlignment', 'right',...
	        'String',e,...
	        'ToolTip','Insert the size of neighbourhood.' );
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h2(3) h1(4)])

  if method==12 | method==9, set(h0,'Enable','Off'); end

  dark_factor=.86;
  h0=uicontrol(props.frame,'BackgroundColor',dark_factor*props.frame.BackgroundColor, ...  % Frame LOSsearch
            'Units','Normalized',...
	        'Position',[.1 .212 .8 .22]);

  h0=uicontrol(props.text,...		% Text LOSsearch
            'Units','Normalized',...
	        'FontAngle','italic', ...
	        'String','LOS Search',...
	        'Tag','TextLOSsearch',...
	        'BackgroundColor',dark_factor*props.frame.BackgroundColor,...
	        'Enable','off',...
	        'Position',[.12 .40 .6 .02]);
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h1(3) h1(4)])

  h0=uicontrol(props.button,...		% Set Point
            'Units','Normalized',...
	        'Tag','SetPoint',...
	        'String','Set Point',...
	        'Position',[.16 .36 .314 .032],...
	        'Enable','off',...
	        'BackgroundColor',dark_factor*props.frame.BackgroundColor,...
	        'ToolTip','Set a Point on LOS.',...
	        'CallBack','crp2 LOSset');

  h0=uicontrol(props.button,...		% Clear Point
            'Units','Normalized',...
	        'Tag','ClearPoint',...
	        'String','Clear P',...
	        'Position',[.52 .36 .314 .032],...
	        'Enable','off',...
	        'ToolTip','Clear a Point on LOS.',...
	        'BackgroundColor',dark_factor*props.frame.BackgroundColor,...
	        'CallBack','crp2 LOSclear');

  h0=uicontrol(props.button,...		% Clear All Point
            'Units','Normalized',...
	        'Tag','ClearAllPoint',...
	        'String','Clear All',...
	        'Position',[.52 .32 .314 .032],...
	        'Enable','off',...
	        'ToolTip','Clear All Points on LOS.',...
	        'BackgroundColor',dark_factor*props.frame.BackgroundColor,...
	        'CallBack','crp2 LOSallclear');

  h0=uicontrol(props.text,...		% Text LOSwidthX
            'Units','Normalized',...
	        'Tag','LOSwidthXtext',...
	        'String','dx:',...
	        'BackgroundColor',dark_factor*props.frame.BackgroundColor,...
	        'Position',[.16 .274 .35 .02],...
	        'Enable','off');
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h1(3) h1(4)])

  h0=uicontrol(props.edit,...		% Input LOSwidthX
            'Units','Normalized',...
	        'Tag','LOSwidthX',...
	        'Position',[.315 .278 .15 .026],...
	        'String','1',...
	        'Enable','off',...
	        'ToolTip','Insert the LOS width in X-direction.' );
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h2(3) h1(4)])


  h0=uicontrol(props.text,...		% Text LOSwidthY
            'Units','Normalized',...
	        'Tag','LOSwidthYtext',...
	        'String','dy:',...
	        'Position',[.52 .274 .35 .02],...
	        'Enable','off',...
	        'BackgroundColor',dark_factor*props.frame.BackgroundColor);
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h1(3) h1(4)])

  h0=uicontrol(props.edit,...		% Input LOSwidthY
            'Units','Normalized',...
	        'Tag','LOSwidthY',...
	        'Position',[.679 .278 .15 .026],...
	        'String','1',...
	        'Enable','off',...
	        'ToolTip','Insert the LOS width in Y-direction.' );
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h2(3) h1(4)])


  h0=uicontrol(props.button,...		% Button ApplyLOSsearch
            'Units','Normalized',...
	        'String','Apply',...
	        'Position',[.16 .23 .314 .032],...
	        'Tag','Apply2',...
	        'Callback','crp2 LOSsearch',...
	        'BackgroundColor',dark_factor*props.frame.BackgroundColor,...
	        'Enable','off',...
	        'ToolTip','Searches the LOS.');


  h0=uicontrol(props.button,...		% Button Stores LOS
            'Units','Normalized',...
	        'String','Store',...
	        'Position',[.52 .23 .314 .032],...
	        'Tag','Store2',...
	        'Callback','crp2 store2',...
	        'BackgroundColor',dark_factor*props.frame.BackgroundColor,...
	        'Enable','off',...
	        'ToolTip','Stores the LOS.');


  if ~isunix
    h0=uicontrol(props.frame, ... % Frame Embedding
            'Units','Normalized',...
	        'Position',[.1 .15 .38 .045]);
  end

  h0=uicontrol(props.checkbox, ...		% Checkbox Stretch Plot
            'Units','Normalized',...
	        'Position',[.1 .15 .38 .045], ...
	        'String','Stretch', ...
	        'Tag','Stretch',...
	        'CallBack','crp2 stretch',...
	        'Value',1,...
	        'ToolTip','Streches the plotted CRP to a squared plot.' );

  h0=uicontrol(props.button, ...		% Button Store Matrix
            'Units','Normalized',...
	        'Position',[.52 .15 .38 .045], ...
	        'String','Store Matrix', ...
	        'Style','pushbutton', ...
	        'Tag','Store',...
	        'CallBack','crp2 store',...
	        'Value',1,...
	        'ToolTip','Stores the CRP matrix into variable X in the workspace.' );

  h0=uicontrol(props.button,...		% Button Help
            'Units','Normalized',...
	        'String','Help',...
	        'Position',[.1 .03 .38 .045],...
	        'Tag','Help',...
	        'Callback','helpwin crp2',...
	        'ToolTip','Opens the helpwindow.');

  h0=uicontrol(props.button,...		% Button Apply
            'Units','Normalized',...
	        'String','Create Recurrence Plot',...
	        'Position',[.1 .09 .8 .045],...
	        'Tag','Apply',...
	        'Callback','crp2 compute',...
	        'ToolTip','Starts the computation - be patient.');

  h0=uicontrol(props.button,...		% Button Close
            'Units','Normalized',...
	        'String','Close',...
	        'Position',[.52 .03 .38 .045],...
	        'Tag','Close',...
	        'Callback','crp2 close',...
	        'ToolTip','Closes CRP windows.');
  set(h9, 'HandleVis','CallBack')

  clear h0 h1 h8 h9

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  vectorexclude

  errcode=4;
  case 'exclude'
	
	ds = get(findobj('Tag','DimEx','Parent',gcf),'UserData');
	j=str2num(get(gco,'String'));
	if ds(j,j)==1
	   set(gco,'ForegroundColor', [.5 0 0],...
	           'FontWeight','bold',...
		   'BackgroundColor', .8*props.button.BackgroundColor,...
		   'ToolTip','Include this vector.');
	   ds(j,j)= -1;
	   if -(trace(ds))==length(ds)
	     ds(j,j)= 1;
	     set(gco,'ForegroundColor', [0 0 0],'FontWeight','normal',...
	           'BackgroundColor', props.button.BackgroundColor,...
		   'ToolTip','Exclude this vector.');
	   else
             h0=gcf; figure(hCRP)
	     h1=findobj('Tag','Data1'); set(h1(length(h1)-j+1),'Visible','off')
	     h2=findobj('Tag','Data2'); set(h2(length(h1)-j+1),'Visible','off')
	     figure(h0), clear h1 h0
	   end
	else
	   set(gco,'ForegroundColor', [0 0 0],'FontWeight','normal',...
	           'BackgroundColor', props.button.BackgroundColor,...
		   'ToolTip','Exclude this vector.');
	   ds(j,j)= 1;
	   h0=gcf; figure(hCRP)
	   h1=findobj('Tag','Data1'); set(h1(length(h1)-j+1),'Visible','on')
	   h2=findobj('Tag','Data2'); set(h2(length(h1)-j+1),'Visible','on')
	   figure(h0), clear h1 h0
	end
	set(findobj('Tag','DimEx','Parent',gcf),'UserData', ds);
	m1=length(find(sum(ds)+1));
	set(findobj('Tag','Dim','Parent',gcf),'string', ['x ',num2str(m1)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% delay

  case 'dimfit'
  
  errcode=5;
  m1=length(find(sum(ds)+1));
  if m0~=1
     set(findobj('Tag','Delay','Parent',gcf),'Enable', 'On');
  else
     set(findobj('Tag','Delay','Parent',gcf),'Enable', 'Off');
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
  h1=findobj('Parent',hCRP,'Tag','CRPPlot');
  h2=findobj('Tag','Diagonal','Parent',h1);
  if ~nogui
      if isempty(h2) 
        h2=line('Parent',h1,'visible','off','Tag','Diagonal','color',[1 0 0],...
        'LineWidth',1);
      else
        set(h2,'visible','off')
      end
      set(findobj('Tag','Store2','Parent',hCtrl),'Enable','Off') 
      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'Visible','on')
      set(findobj('Tag','CRPData','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'Visible','off')
  end

  % check if plugin exist and is executable
  [plugin_exist, plugin_name, plugin_path] = is_crp_plugin;
  if plugin_exist & ( mflag < 4 | mflag == 9 | mflag == 10 ) & length(x) == length(y) & ~ispc % if plugin exist and method is MAX, MIN, EUC ord DIS
      
      if nogui == 1, disp('(plugin used)'), end
      [X matext] = crp_plugin(x, y, m0, t, e, mflag, hCRP, plugin_path, 0);

  else
  % use builtin implementation
      if ~nogui
          set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Import Embedding Vectors'),drawnow
      end
      ex=ceil(find(~(ds+1))/m);
      x(:,ex)=[]; y(:,ex)=[]; 
      m=m-length(ex);
      if m==0, 
        set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','No embedding vectors!'),drawnow
        return
      end

      if m0>1
        x2=x(1:end-t*(m0-1),:);
        y2=y(1:end-t*(m0-1),:);
        for i=1:m0-1,
          if check_stop(hCRP,hCtrl,nogui,obj), break, end
          x2(:,m*i+1:m*(i+1))=x(1+t*i:end-t*(m0-i-1),:);
          y2(:,m*i+1:m*(i+1))=y(1+t*i:end-t*(m0-i-1),:);
        end
        x=x2; y=y2; Nx=size(x,1); Ny=size(y,1);
        m=m0*m; clear x2 y2
      end

      x1=repmat(x,1,Ny);
      if check_stop(hCRP,hCtrl,nogui,obj), return, end
      for mi=1:m, x2(:,mi)=reshape(rot90(x1(:,0+mi:m:Ny*m+mi-m)),Nx*Ny,1); end
      y1=repmat(y,Nx,1); x1=x2; clear x2
      if check_stop(hCRP,hCtrl,nogui,obj), return, end

      if ~nonorm, unit = ''; else unit = '\sigma'; end

      switch(mflag)


      %%%%%%%%%%%%%%%%% local CRP, fixed distance maximum norm

      case {1,2,3,5}

      errcode=111;
      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Compute Distance Matrix'),drawnow
      s=zeros(m,Nx*Ny);
      if check_stop(hCRP,hCtrl,nogui,obj), return, end
      s1=zeros(Nx*Ny,m);
      if check_stop(hCRP,hCtrl,nogui,obj), return, end
      X=uint8(zeros(Ny,Nx));
      if check_stop(hCRP,hCtrl,nogui,obj), return, end
      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Compute Distance Matrix'),drawnow
      s1=x1-y1;
      switch mflag
        case 1
      %%%%%%%%%%%%%%%%% maximum norm
        if m>1
          if check_stop(hCRP,hCtrl,nogui,obj), return, end
          s=max(abs(s1'));
        else
          s=abs(s1);
        end
        matext=[num2str(round(100*e)/100) unit ' (fixed distance maximum norm)'];
        case 5
      %%%%%%%%%%%%%%%%% maximum norm, fixed RR
          errcode=115;
          if m>1
            if check_stop(hCRP,hCtrl,nogui,obj), return, end
            s=max(abs(s1'));
          else
            s=abs(s1);
          end
          ss = sort(s(:));
          idx = ceil(e * length(ss));
          e = ss(idx);
          matext=[num2str(round(100*e)/100) '\sigma (fixed distance maximum norm, fixed RR)'];
        case 2
      %%%%%%%%%%%%%%%%% euclidean norm
        errcode=112;
        if m>1
          if check_stop(hCRP,hCtrl,nogui,obj), return, end
          set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Transpose  Matrix'),drawnow
          s=s1.^2';
          if check_stop(hCRP,hCtrl,nogui,obj), return, end
          set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Sum Matrix'),drawnow
          s1=sum(s);
          if check_stop(hCRP,hCtrl,nogui,obj), return, end
          set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Squareroot Matrix'),drawnow
          s=sqrt(s1);
        else
          s=abs(s1);
        end
        matext=[num2str(round(100*e)/100) unit ' (fixed distance euclidean norm)'];
        case 3
      %%%%%%%%%%%%%%%%% minimum norm
        errcode=113;
        if m>1
          if check_stop(hCRP,hCtrl,nogui,obj), return, end
          s = sum(abs(s1),2);
        else
          s=abs(s1);
        end
        matext=[num2str(round(100*e)/100) unit ' (fixed distance minimum norm)'];
      end
      if check_stop(hCRP,hCtrl,nogui,obj), return, end

      %%%%%%%%%%%%%%%%% now apply threshold
      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')'),'String','Building CRP Matrix'),drawnow
      X2 = s<e;
      X=reshape(uint8(X2),Ny,Nx); 
      clear s s1 x1 y1



      %%%%%%%%%%%%%%%%% local CRP, normalized distance euclidean norm

      case 4

      errcode=114;
      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Normalize Embedding Vectors'),drawnow
      Dx=sqrt(sum(((x1(:,:)).^2)'))';
      if check_stop(hCRP,hCtrl,nogui,obj), return, end
      Dy=sqrt(sum(((y1(:,:)).^2)'))';
      if check_stop(hCRP,hCtrl,nogui,obj), return, end
      x1=x1./repmat(Dx,1,m);
      if check_stop(hCRP,hCtrl,nogui,obj), return, end
      y1=y1./repmat(Dy,1,m); clear Dx Dy 
      if check_stop(hCRP,hCtrl,nogui,obj), return, end

      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Reshape Embedding Vectors'),drawnow

      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Compute Distance Matrix'),drawnow
      if m>1
         s=sqrt(sum(((x1(:,:)-y1(:,:)).^2)'));
      else
         s=abs(x1(:,1)-y1(:,1));
      end
      if check_stop(hCRP,hCtrl,nogui,obj), return, end

      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Building CRP Matrix'),drawnow
      X=uint8(reshape((s/max(s))<(e/max(s)),Ny,Nx)); clear s x1 y1 
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

      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Reshape Embedding Vectors'),drawnow

      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Compute Distance Matrix'),drawnow
      if m>1
         s=sqrt(sum(((x1(:,:)-y1(:,:)).^2)'));
      else
         s=abs(x1(:,1)-y1(:,1));
      end
      if check_stop(hCRP,hCtrl,nogui,obj), return, end

      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Sort Distance Matrix'),drawnow
      mine=round(Ny*e);
      [SS, JJ]=sort(reshape(s,Ny,Nx)); JJ=JJ';
      if check_stop(hCRP,hCtrl,nogui,obj), return, end

      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Building CRP Matrix'),drawnow
      X1(Nx*Ny)=uint8(0); X1(JJ(:,1:mine)+repmat([0:Ny:Nx*Ny-1]',1,mine))=uint8(1);
      if check_stop(hCRP,hCtrl,nogui,obj), return, end
      X=reshape(X1,Ny,Nx); clear X1 SS JJ s

      matext=[num2str(round(1000*mine/Ny)/10) '% (fixed neighbours amount)'];


      %%%%%%%%%%%%%%%%% local CRP, interdependent neigbours 

      case 7

      errcode=117;
      warning('Method not available!')
      matext='';


      %%%%%%%%%%%%%%%%% order matrix

      case 8

      errcode=118;
      warning('Method not available!')
      matext='';


      %%%%%%%%%%%%%%%%% order pattern recurrence plot

      case 9

      errcode=119;

      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Compute Order Patterns'),drawnow

      % create a order pattern test
      cmdStr = '';
      for i=2:m;
          cmdStr = [cmdStr, ' permX(:,', num2str(i-1) ,') < permX(:,', num2str(i), ') + eps'];
          if i < m, cmdStr = [cmdStr, ' &']; end
      end
      if m==1
          cmdStr = '1'; 
          disp('Warning: No order patterns for dimension one; please use higher dimension!')
      end

      % order patterns by permutation of the set of values
      clear patt*
      for i=1:size(x,1)
          permX=perms(x(i,:));
          orderPattern = find(eval(cmdStr));
          if isempty(orderPattern) orderPattern = 0; end
          pattX(i) = orderPattern(1);
      end
      for i=1:size(y,1)
          permX=perms(y(i,:));
          orderPattern = find(eval(cmdStr));
          if isempty(orderPattern) orderPattern = 0; end
          pattY(i) = orderPattern(1);
      end

      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Create Order Patterns Matrix'),drawnow

      px = permute(pattX', [ 1 3 2 ]);
      py = permute(pattY', [ 3 1 2 ]);
      X = uint8(px(:,ones(1,length(y)),:) == py(ones(1,length(x)),:,:));
      X = X';
      if check_stop(hCRP,hCtrl,nogui,obj), return, end

      matext='';

      %%%%%%%%%%%%%%%%% Levenshtein

      case 10

      errcode=120;

      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Compute Distance Matrix'),drawnow
      s1 = levenshtein(x1, y1);
      X = reshape(s1,Ny,Nx);
      matext='';


      %%%%%%%%%%%%%%%%% DTW

      case 11

      errcode=121;

      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Compute Distance Matrix'),drawnow
      s1 = dtw(x1, y1);
      X = reshape(s1,Ny,Nx);
      matext='';



      %%%%%%%%%%%%%%%%% global CRP

      case length(check_meth)

      errcode=110+length(check_meth);
      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Reshape Embedding Vectors'),drawnow

      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Compute Distance Matrix'),drawnow
      if m>1
         s=sqrt(sum(((x1(:,:)-y1(:,:)).^2)'));
      else
         s=abs(x1(:,1)-y1(:,1));
      end
      if check_stop(hCRP,hCtrl,nogui,obj), return, end

      set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Building CRP Matrix'),drawnow
      X=reshape(s,Ny,Nx);
      matext='';

      end

  end % end plugin

  if nogui==0
    for i=1:length(obj.enable), set(obj.children(i),'Enable',obj.enable{i}); end
    set(h(1),'String','Create Recurrence Plot',...
             'ToolTip','Starts the computation - be patient.',...
             'Callback','crp2 compute')
    setptr([hCRP,hCtrl],'arrow')
  end

  %%%%%%%%%%%%%%%%% show CRP 

  if nogui==0
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
   if nargout==1, xout=X; end
 end

 if strcmpi(get(findobj('Tag','SetPoint','Parent',hCtrl),'Enable'),'off')
     set(findobj('Tag','TextLOSsearch','Parent',hCtrl),'Enable','On')
     set(findobj('Tag','SetPoint','Parent',hCtrl),'Enable','On')
     set(findobj('Tag','ClearPoint','Parent',hCtrl),'Enable','On')
     set(findobj('Tag','ClearAllPoint','Parent',hCtrl),'Enable','On')
     set(findobj('Tag','LOSwidthXtext','Parent',hCtrl),'Enable','On')
     set(findobj('Tag','LOSwidthX','Parent',hCtrl),'Enable','On')
     set(findobj('Tag','LOSwidthYtext','Parent',hCtrl),'Enable','On')
     set(findobj('Tag','LOSwidthY','Parent',hCtrl),'Enable','On')
     set(findobj('Tag','Apply2','Parent',hCtrl),'Enable','On')
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% store2

  case 'store2'
  
  errcode=15;
  t=get(findobj('Tag','Apply2','Parent',hCtrl),'UserData');
  if isempty(t)
     warndlg('The LOS vector is still empty. Please start the computation of the LOS before storing.','No LOS')
     waitforbuttonpress
  else
     assignin('base','t_out', t)
  end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOSmove

  case 'LOSmove'
  
  errcode=16;
  if isempty(get(gco,'UserData'))
    set(gco,'UserData',1,'ButtonDownFcn','crp2 LOSmove_end')
    set(gcf,'WindowButtonMotionFcn','crp2 LOSmove')
  end
  h1 = round(get(gca,'CurrentPoint'));
  set(gco, 'XData', h1(1,1), 'YData', h1(1,2))
  clear h1

  case 'LOSmove_end'
  errcode=161;
  if ~isempty(get(gco,'UserData'))
    set(gco,'UserData',[],'ButtonDownFcn','')
    set(gcf,'WindowButtonMotionFcn','')
  end
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOSclear

  case 'LOSclear'
  
  errcode=17;
  if gcf==hCtrl
    figure(hCRP)
    k=waitforbuttonpress;
    h1 = round(get(gca,'CurrentPoint')); h1(2,:)=[]; h1(:,3)=[];
    finalRect = rbbox;
    h2 = round(get(gca,'CurrentPoint')); h2(2,:)=[]; h2(:,3)=[];
    h(1,1)=min(h1(:,1),h2(:,1));h(2,1)=max(h1(:,1),h2(:,1));
    h(1,2)=min(h1(:,2),h2(:,2));h(2,2)=max(h1(:,2),h2(:,2));
    
    i=find(fixp(:,1)>=h(1,1) & fixp(:,2)>=h(1,2) & fixp(:,1)<=h(2,1) & fixp(:,2)<=h(2,2));
    for j=1:length(i); delete(findobj('tag','fixp','Parent',findobj('Parent',hCRP,'Tag','CRPPlot'),'xdata',fixp(i(j),1))), end
  else
    h(1,1)=get(gco,'XData');h(1,2)=get(gco,'YData');
    delete(gco)
  end
  clear h

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOSallclear

  case 'LOSallclear'
  
    errcode=171;
    delete(findobj('tag','fixp','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')))
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOSset

  case 'LOSset'
  
  errcode=18;
  if gcf==hCtrl
    figure(hCRP)
    ginput(1);
  end
  
  h=round(get(gca,'currentp'));
  fixp(end+1,1)=h(1,1); fixp(end,2)=h(1,2);
  [i j]=sort(fixp(:,1));
  fixp=fixp(j,:);
  h0=uicontextmenu;
  line(h(1,1),h(1,2),1000,...
       'MarkerSize',12,...
       'Marker','.',...
       'Color',[1 0 0],...
       'Tag','fixp',...
       'UIContextMenu',h0)
  uimenu(h0, 'Label', 'Set Point', 'Callback', 'crp2 LOSset')
  uimenu(h0, 'Label', 'Move Point', 'Callback', 'crp2 LOSmove')
  uimenu(h0, 'Label', 'Clear Point', 'Callback', 'crp2 LOSclear')
  clear h h0

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOSsearch

  case 'LOSsearch'

  errcode=19;
  if nogui==0
    setptr([hCRP,hCtrl],'watch')
    h=findobj('Tag','Apply2','Parent',hCtrl);
    obj.children=get(hCtrl,'Children');
    obj.enable=get(obj.children,'Enable'); 
    set(obj.children,'Enable','off')
    set(h(1),'String','Stop',...
             'ToolTip','Stops the computation.',...
             'Enable','on',...
	     'Callback','set(gcbo,''String'',''Stopped'')')
  end
  x0=0;
  y0=0;
  Dmax1=str2num(get(findobj('Tag','LOSwidthX','Parent',hCtrl),'string'));
  Dmax2=str2num(get(findobj('Tag','LOSwidthY','Parent',hCtrl),'string'));
  X=double(get(findobj('Tag','CRPData','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'UserData'));
  
  N=size(X); Nleer=0; Nvoll=0;

  minvalue=min(min(X));maxvalue=max(max(X));
  if (length(find(X==minvalue))+length(find(X==maxvalue))==length(X(:)))

  flagpoint=0;
   
% looks for the beginning of the diagonal LOS

  errcode=191;
  for i=1:N(2);
    if check_stop_LOS(hCRP,hCtrl,nogui,obj), break, end
    if i<=N(1)
      if ~isempty(fixp), if i>=fixp(1,1), x0=fixp(1,1); y0=fixp(1,2); end, end
      if X(i+y0,1+x0)~=0 
        v=i+y0;
        h=1+x0;
        break
      end
    end
    if X(1+y0,i+x0)~=0
      h=i+x0;
      v=1+y0;
      break
    end
  end
  
  Nleer=i-1;
%  t(1:v)=h;
  t(1+x0:h)=v;

% start estimation of the LOS

  errcode=192;
  mflag=0;
  i=2;
%  h_wait=waitbar(0,'Compute LOS ...',props.window);
  h_wait=waitbar(0,'Compute LOS ...');
  set(h_wait,'HandleVisibility','on');

  while h<N(2)-1 & v<N(1)-1 & mflag~=1,
    waitbar(h/N(2))
    if check_stop_LOS(hCRP,hCtrl,nogui,obj), break, end

    dw=1;

    if v < 1, v = 1; end
    if h < 1, h = 1; end
  
    W=X(v:v+dw,h:h+dw);
    W(1,1)=0; 

    if ~isempty(fixp)
      if find(h==fixp(:,1))
        dv=fixp(find(h==fixp(:,1)),2)-v;
        dh=1; 
        if dv <= 0
            h0 = min(find(t > v+dv));
            if isempty(h0), h0 = h; end
            t0 = t(end) + dv;
            fixold = find(h==fixp(:,1));
            if ~isempty(fixold) & fixold > 1; fixold = fixp(fixold - 1,1); else fixold = 0; end
            h1 = max([h0-min([Dmax1,Dmax2]);1;fixold]);
            h2 = max([h0;1;fixold]);
%            h1 = max([h-2*(h-h0);1])
            if length(t) > h1 & length(t) > h2
                t1 = t(h1);
                t2 = t(h2);
                p = polyfit([h1 h2 h],[t1 t2 t0],1);
                t(h1:h) = polyval(p,h1:h);
%                t(h1:h) = spline([h1 h2 h],[t1 t2 t0],h1:h);
            end
            v = v + dv;
            dv = 1; dw = 0;
        end
        flagpoint=1;
      end
    end

    if flagpoint==0;

% looks for the existence of the next recurrence point in diagonal direction

    errcode=193;
    while sum(sum(W))==0 & mflag==0 & flagpoint==0,
        Nleer=Nleer+1;
        if ~isempty(fixp)
           if find(h+dw==fixp(:,1))
  	           W=1; dh=dw;
	           flagpoint=1;
	           break
           end
        end
        dw=dw+1;
        if v+dw < N(1) & h+dw < N(2)
           W=X(v:v+dw,h:h+dw);
           W(1,1)=0; 
        else
           mflag=1;
        end
    end

    if mflag==1
      break
    end

% determines the coordinates of the next recurrence point
  
    errcode=194;
    dh0=min(find(sum(W)));
    dv0=min(find(sum(W')));
    dh=min(find(W(dv0,:)));
    dv=min(find(W(:,dh0)));
    if dh>dh0, dh=dh0;end
    if dv>dv0, dv=dv0;end
  

% determines the local width of the diagonal LOS

    errcode=195;
    dh1=dh;
    dv1=dv;
    
    
    % neues Fenster
    WL1=Dmax1;
    WL2=Dmax2;
    if v+dv1-2+WL2>=N(1), WL2=N(1)-(v+dv1-2); end
    if h+dh1-2+WL1>=N(2), WL1=N(2)-(h+dh1-2); end
    Wn=X(v+dv1-1:v+dv1-2+WL2,h+dh-1:h+dh-2+WL1);
    % Schwerpunkt davon ausrechnen
    if sum(sum(Wn))~=0, Sh=sum(sum(Wn).*[1:WL1])/sum(sum(Wn));else Sh=WL1/2;end
    if sum(sum(Wn'))~=0,Sv=sum(sum(Wn').*[1:WL2])/sum(sum(Wn'));else Sv=WL2/2;end
    % neue Dmax berechnen
    if Sh>=Sv, Dmax2n=WL1*Sv/Sh; Dmax1n=WL1; end
    if Sv>Sh,  Dmax1n=WL2*Sh/Sv; Dmax2n=WL2; end

%    while X(v+dv1-1,h:dh-1)==1 & v+dv1-1<N(1) &dv1<Dmax2

    while X(v+dv1-1,h+dh-1)==1 & v+dv1-1<N(1) &dv1<Dmax2n
     dv1=dv1+1;
    end
%    while X(v:dv-1,h+dh1-1)==1 & h+dh1-1<N(2) & dh1<Dmax1
    while X(v+dv-1,h+dh1-1)==1 & h+dh1-1<N(2) & dh1<Dmax1n
      if ~isempty(fixp)
        if find(h+dh1-1==fixp(:,1)) & find(h+fix((dh1+dh)/2)==fixp(:,1))
          dv=fixp(find((h+dh1-1)==fixp(:,1)),2)-v;
	  if isempty(dv), dv=fixp(find((h+fix((dh1+dh)/2))==fixp(:,1)),2)-v; end
          flagpoint=1;
          break
        end
      end
      dh1=dh1+1;
    end


% compute the mean of the diagonal LOS
  
    errcode=196;
    if flagpoint==0
      dh=fix((dh1+dh)/2);
      dv=fix((dv1+dv)/2);

      if dh>0
        dh=dh-1;
      end
      if dv>0
        dv=dv-1;
      end
%      dh=dh+dh1;
    end
    
    end % flagpoint end
    flagpoint=0;


% output
    if dh~=0 & dv~=0
      t(h:h+dh)=v:dv/dh:v+dv;
    elseif dv~=0
      t(h)=v+dv;
    elseif dv==0
      t(h:h+dh)=v;
    end
    Nvoll=Nvoll+sqrt(dv^2+dh^2);

    if dh==0 & dv==0
      dh=1;
    end

% moves the startpoint for further looking
    h=h+dh;
    v=v+dv;

  end

  else
% DTW algorithm

   h_wait=waitbar(0,'Compute LOS ...');
   set(h_wait,'HandleVisibility','on',props.window);
  
   t=1;
   i=y0+1; j=x0+1;
   flag2=0;
   while i<N(1)-Dmax1-1 & j<N(2)-Dmax2-1
     errcode=197;
     waitbar(i/N(1)), j0=j;
%     [temp pos]=min([sum(sum(X(i:i+Dmax1,j+1:j+1+Dmax2))),sum(sum(X(i+1:i+1+Dmax1,j+1:j+1+Dmax2))),sum(sum(X(i+1:i+1+Dmax1,j:j+Dmax2)))]);
     [temp pos]=min([(mean(X(i,j+1:j+1+Dmax2))),(mean(diag(X(i+1:i+1+Dmax1,j+1:j+1+Dmax2)))),(mean(X(i+1:i+1+Dmax1,j)))]);
     switch(pos)
       case 1
        j=j+1;
	flag1=1;
       case 2
        i=i+1; j=j+1;
	flag1=1;
	flag2=1;
       case 3
        i=i+1;
	flag2=1;
     end
     if ~isempty(fixp) & flag1==1 & flag2==1
       errcode=1981;
       h=find(fixp(:,1)==j);
       if ~isempty(h) & i<max(fixp(h,2)); 
       h2=find(i<=fixp(h,2)); 
       i=fixp(h(h2(1)),2);  flag1=0; end
     end
     if ~isempty(fixp) & flag2==1 & flag1==1
       errcode=1982;
       h=find(fixp(:,2)==i);
       if ~isempty(h) & j<max(fixp(h,1)); 
       h2=find(j<=fixp(h,1));
       j=fixp(h(h2(1)),1);
       t(j0:j-1)=t(j0);
       flag2=0; end
     end
     t(j)=i;

   end

   t(find(~t))=1;

  end

  delete(h_wait)
  errcode=199;
  
  h1=findobj('Parent',hCRP,'Tag','CRPPlot');
  h2=findobj('Tag','Diagonal','Parent',h1);
  if isempty(h2)
    h2=line('Parent',h1,'visible','on','Tag','Diagonal','color',[1 0 0],...
    'LineWidth',1);
  end
  set(h2,'xdata',[1:length(t)],'ydata',t,'visible','on')
  set(findobj('Tag','Apply2','Parent',hCtrl),'UserData',t)
  if nogui==0
    setptr([hCRP,hCtrl],'arrow')
    h=findobj('Tag','Apply2','Parent',hCtrl);
    set(h(1),'String','Apply',...
     	         'ToolTip','Searches the LOS.',...
	         'Callback','crp2 LOSsearch')
    for i=1:length(obj.enable), set(obj.children(i),'Enable',obj.enable{i}); end
  end
  set(findobj('Tag','Store2','Parent',hCtrl),'Enable','On') 

end
try, set(0,props.root), end
warning on
set(0,'ShowHidden','Off')

%%%%%%% error handling

%if 0
catch
  try
    if nogui==0
        for i=1:length(obj.enable), set(obj.children(i),'Enable',obj.enable{i}); end
        set(h(1),'String','Apply',...
                 'ToolTip','Starts the computation - be patient.',...
                 'Callback','crp compute')
        setptr([hCRP,hCtrl],'arrow')
    end
    if nargout, xout = NaN; end
  end
  z=whos;x=lasterr;y=lastwarn;in=varargin{1};
  print_error('crp2',z,x,y,in,mflag,action)
  try, set(0,props.root), end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=check_stop_LOS(hCRP,hCtrl,nogui,obj)

    global errcode
    errcode=errcode+.02;
    out=0;
    if nogui==0
      h=findobj('Tag','Apply2','Parent',hCtrl);
      if strcmpi(get(h(1),'String'),'stopped')
        set(h(1),'String','Apply',...
     	         'ToolTip','Searches the LOS.',...
	         'Callback','crp2 LOSsearch')
        for i=1:length(obj.enable), set(obj.children(i),'Enable',obj.enable{i}); end
        set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')'),'String','Stopped'),drawnow
        setptr([hCRP,hCtrl],'arrow')
        out=1;
      end
    end
